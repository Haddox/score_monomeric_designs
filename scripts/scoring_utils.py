"""
Extra modules for scoring protein structures
"""

import os
import pandas
import doctest
import subprocess
import itertools
import time
import pyrosetta
import pyrosetta.rosetta
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from pyrosetta.rosetta.core.select.residue_selector import LayerSelector
from pyrosetta.bindings.utility import bind_method

import Bio.PDB
BIO_PDB_parser = Bio.PDB.PDBParser(QUIET=True)

@bind_method(HBondSet)
def __iter__(self):
    for i in range(1, self.nhbonds() + 1):
        yield self.hbond(i)

@bind_method(HBondSet)
def hbonds(self):
    return [hb for hb in self]

# Set up a score function in `PyRosetta`
pyrosetta.init('-beta_nov16 -corrections:beta_nov16 -mute all -mute core -mute protocols')
sf = pyrosetta.get_fa_scorefxn()
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(sf)


# Define modules
def ComputeAminoAcidFrequencyInResidueSubset(aa_seq, per_residue_ids, aa, subset_id, return_counts=False):
    """
    Compute the frequency of a given amino acid among a subset of residues with a common ID

    Args:
        `aa_seq` (string) : the amino-acid sequence with single-letter amino-acid
            identifiers
        `per_residue_ids` (string) : a string of characters that is the same length
            as `aa_seq`, providing an ID to each of the corresponding amino acids.
            For instance, the string 'HEH' would indicate that the first and third
            residues have an ID equal to 'H', while the second residue has an ID
            equal to 'E'.
        `aa` (string) : the function will compute the frequency of this amino acid
        `subset_id` (string) : only the subset of residues with this ID will be
            considered when computing the frequency of the given amino acid
        `return_counts` (bool, False) : if True, this function will return the
            counts of an amino acid instead of its frequency

    Returns:
        By default, this function returns the frequency (float) of the specified amino
            acid among all residues that have a common ID equal to `subset_id`.
            However, if `return_counts` is changed to `True`, then it will return the
            counts (int) of the amino acid instead of the frequency. If the `subset_id`
            does not occur in the `aa_seq`, it will return `None`

    Code for `doctest`:
        >>> aa_seq = 'ACDAAAC'
        >>> per_res_ids = 'HHHEEEE'
        >>> freq = ComputeAminoAcidFrequencyInResidueSubset(aa_seq, per_res_ids, 'A', 'E')
        >>> round(float(freq), 2)
        0.75
        >>> ComputeAminoAcidFrequencyInResidueSubset(aa_seq, per_res_ids, 'A', 'E', return_counts=True)
        3
        >>> freq = ComputeAminoAcidFrequencyInResidueSubset(aa_seq, per_res_ids, 'D', 'H')
        >>> round(float(freq), 2)
        0.33
        >>> freq = ComputeAminoAcidFrequencyInResidueSubset(aa_seq, per_res_ids, 'D', 'L')
        >>> round(float(freq), 2)
        0.0
    """

    # Conver the strings into pandas series
    assert len(aa_seq) == len(per_residue_ids), \
        "The amino-acid sequence and per-res secondary structures must have the same lengths"
    aa_seq = pandas.Series(list(aa_seq))
    per_residue_ids = pandas.Series(list(per_residue_ids))

    # Identify which residues have the specified ID, returning `None` if there
    # aren't any residues with the ID
    residues_in_subset_bools = per_residue_ids == subset_id
    n_residues_in_subset = sum(residues_in_subset_bools)
    if n_residues_in_subset == 0:
        return None
    else:
        # Compute the counts or frequency of a given amino acid among the subset
        # of residues defined above and return the result
        aa_count_in_subset = aa_seq[
            (residues_in_subset_bools) & \
            (aa_seq == aa)
        ].count()
        if return_counts:
            return aa_count_in_subset
        else:
            return aa_count_in_subset / float(n_residues_in_subset)


def find_hbonds(p, derivatives=False, exclude_bb=True, exclude_bsc=True,
                exclude_scb=True, exclude_sc=False):
    """Find all hydrogen bonds of a particular type in a supplied Pose
    and return them as a HBondSet.

    Args:
        p (pyrosetta.Pose): The Pose from which to extract hydrogen
            bonds.
        derivatives (bool, optional): Evaluate energy derivatives and
            store them in the HBondSet. Defaults to False.
        exclude_bb (bool, optional): If True, do not store
            backbone--backbone hydrogen bonds in the HBondSet. Defaults
            to True.
        exclude_bsc (bool, optional): If True, do not store
            backbone--side chain hydrogen bonds in the HBondSet.
            Defaults to True.
        exclude_scb (bool, optional): If True, do not store
            side chain--backbone hydrogen bonds in the HBondSet.
            Defaults to True.
        exclude_sc (bool, optional): If True, do not store
            side chain--side chain hydrogen bonds in the HBondSet.
            Defaults to False.

    Returns:
        pyrosetta.rosetta.core.scoring.hbonds.HBondSet: A hydrogen bond
        set containing the specified types of hydrogen bonds in the
        Pose.
    """
    p.update_residue_neighbors()
    hbset = HBondSet()
    fill_hbond_set(
        p, derivatives, hbset, exclude_bb, exclude_bsc, exclude_scb,
        exclude_sc
    )
    return hbset


def layer_selector_mover(layer, core_cutoff=3.7, surface_cutoff=1.3):
    """
    Set up a PyRosetta Mover that can be used to select a specific layer using the
    side-chain neighbor algorithm

    Args:
        `layer` (string) : the layer to be selected. This variable can be "core",
            "boundary", or "surface".
        `core_cutoff` (float) : the cutoff used to define the core using the
            side-chain neighbor algorithm. Residues with at least this many
            neighbors are considered to be in the core.
        `surface_cutoff` (float) : the cutoff used to define the surface using
            the side-chain neighbor algorithm. Residues with fewer than this
            many neighbors are considered to be on the surface.

    Returns:
        `select_layer` (PyRosetta mover) : a PyRosetta LayerSelector Mover that
            can be applied be applied to a pose to select the layer specified by
            the input arguments.
    """
    # Initiate the mover, defining layer cutoffs
    select_layer = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    select_layer.set_use_sc_neighbors(True)
    select_layer.set_cutoffs(core=core_cutoff, surf=surface_cutoff)

    # Select the specific layer given the input
    layer_bools = {
        'core' : [True, False, False],
        'boundary' : [False, True, False],
        'surface' : [False, False, True]
    }
    pick_core, pick_boundary, pick_surface = layer_bools[layer]
    select_layer.set_layers(
        pick_core=pick_core, pick_boundary=pick_boundary, pick_surface=pick_surface
    )
    return select_layer


def identify_bonded_pairs(p, hbset, layer_selector_mover):
    """Construct an HBondSet that only contains hydrogen bonded
    residues that are within the indicated layer of the supplied
    Pose and return it.

    Args:
        p (pyrosetta.Pose): The Pose to be examined.
        hbset (pyrosetta.rosetta.core.scoring.hbonds.HBondSet): A
            hydrogen bond set to be subsetted based on the burial of
            each residue.
        `layer_selector_mover` (PyRosetta mover) : A mover returned
            by the above function called `LayerSelectorMover`

    Returns:
        pyrosetta.rosetta.core.scoring.hbonds.HBondSet: A hydrogen bond
        set containing only residues on the indicated layer of the Pose.
    """
    return HBondSet(hbset, layer_selector_mover.apply(p))


def find_buried_hbonds(pose, vsasa_burial_cutoff=0.1):
    """
    Find buried hydrogen bonds in an input pose using VSASA to determine burial

    Args:
        `pose`: a `PyRosetta` pose
        `vsasa_burial_cutoff`: the VSASA cutoff below which a residue
            is considered buried. The `BuriedUnsatHbonds` also uses VSASA
            with the same default (default: 0.1).

    Retruns:
        If the pose does not have any buried hydrogen bonds, this function
        returns a value of `None`. But, if it does, this function returns
        a list of tupples, with one tupple per buried hydrogen bond. Each
        tupple has the folling five features:
            `acc_info`: information on the acceptor atom
            `don_info`: information on the donor atom
            `acc_vsasa`: VSASA of the acceptor atom
            `don_vsasa`: VSASA of the donor hydrogen atom
            `don_heavy_vsasa`: VSASA of the donor heavy atom
    """

    # Find all hydrogen bonds in a given pose
    hbset = find_hbonds(
        p=pose, derivatives=False, exclude_bb=True, exclude_bsc=True,
        exclude_scb=True, exclude_sc=False
    )

    # Compute the VSASA of every atom in the pose
    sasa_calc = \
        pyrosetta.rosetta.protocols.vardist_solaccess.VarSolDistSasaCalculator()
    sasa_map = sasa_calc.calculate(pose)

    # Cycle through each hydrogen bond in a pose and identify
    # buried hydrogen bonds
    buried_hbonds = []
    for hb in hbset.hbonds():

        # Get info on the acceptor residue/atom
        acc_res_n = hb.acc_res()
        acc_atom_n = hb.acc_atm()
        acc_res = pose.residues[acc_res_n]
        acc_res_name = acc_res.name3()
        acc_atom_name = acc_res.atom_name(acc_atom_n)
        acc_vsasa = sasa_map[acc_res_n][acc_atom_n]

        # ... and for the donor residue/atom
        don_res_n = hb.don_res()
        don_atom_n = hb.don_hatm()
        don_res = pose.residues[don_res_n]
        don_res_name = don_res.name3()
        don_atom_name = don_res.atom_name(don_atom_n)
        don_vsasa = sasa_map[don_res_n][don_atom_n]

        # ... and for the donor heavy atom that is attached to the polar
        # hydrogen involved in the h bond. I identify the donor heavy atom
        # using the `.bonded_neighbor()` method
        neighbor_atoms = don_res.bonded_neighbor(don_atom_n)
        assert len(neighbor_atoms) == 1
        don_heavy_atom_n = list(neighbor_atoms)[0]
        don_heavy_vsasa = sasa_map[don_res_n][don_heavy_atom_n]

        # Print info about the H bond if it is buried
        if (acc_vsasa < vsasa_burial_cutoff) and \
            (don_heavy_vsasa < vsasa_burial_cutoff):

            acc_info = acc_res_name + '/' + str(acc_res_n) + '/' + acc_atom_name
            don_info = don_res_name + '/' + str(don_res_n) + '/' + don_atom_name
            buried_hbonds.append(','.join(list(map(str, [
                acc_info,
                don_info,
                acc_vsasa,
                don_vsasa,
                don_heavy_vsasa
            ]))))
    if buried_hbonds:
        return buried_hbonds
    else:
        return None

def get_residue_neighbors(pdb_file_name, res_n, distance_cutoff):
    """
    Get a list of residues making side-chain contacts with an
    input residue.

    Specifically, this function iterates over all side-chain
    atoms in the target residue and finds all neighboring
    residues with at least one side-chain atom that is within
    a given distance cutoff of the target residue.

    Args:
        `pdb_file_name`: a path to the input PDB file
        `res_n`: the number of the residue to be used as the
            focal point when selecting surrounding neighbors
        `distance_cutoff`: the distance cutoff to use when
            selecting neighbors.
    Returns:
        A list of residues neighboring the input residue, where
            each residue is listed by its number, and where this
            list includes the number of the input residue.
    """

    # Read in the structure, specify the target residue, and
    # make an object with all atoms in the structure
    structure = BIO_PDB_parser.get_structure('pdb', pdb_file_name)
    assert len(structure) == 1, len(structure)
    chains = Bio.PDB.Selection.unfold_entities(structure, 'C')
    assert len(chains) == 1
    target_residue = structure[0][chains[0].get_id()][res_n]
    atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')

    # Iterate over all atoms in the residue side chain and
    # get a list of all residues that come within a given
    # distance cutoff of each atom
    neighboring_residues = []
    bb_atom_ids = [
        'N', 'H', 'CA', 'HA', 'C', 'O', '1H', '2H', '3H'
    ]
    for target_atom in target_residue.get_atoms():

        # Continue if the atom is a backbone atom
        if target_atom.get_id() in bb_atom_ids:
            continue

        # Make a `NeighborSearch` class object with all atoms in
        # the structure for the query. Then, search for residues
        # with atoms within a given distance cutoff of the target
        # atom, and append them to a growing list
        all_atom_neighbor_search_class = Bio.PDB.NeighborSearch(atoms)
        neighboring_atoms = all_atom_neighbor_search_class.search(
            center=target_atom.coord,
            radius=distance_cutoff,
            level='A'
        )
        neighboring_residues.extend([
            atom.get_parent().get_id()[1]
            for atom in neighboring_atoms
            if atom.get_id() not in bb_atom_ids
        ])

    # Return a list of unique residue numbers of all the
    # contacting residues
    return list(set(neighboring_residues))


def score_design_from_command_line(
    rosettapath, pdb, xml, nstruct, output_dir, filter_name
):
    """
    Score an input PDB using an input XML, then report the results

    Args:
        *rosettapath*: the path to a `rosetta_scripts` app
        *pdb*: the path to the input PDB
        *xml*: the path to the input RosettaScripts XML
        *nstruct*: the number of times to run the protocol
        *output_dir*: the output directory where results will
            be stored
        *filter_name*: the filter to report

    Returns:
        The mean score of the specified filter over all runs
    """

    # Score the PDB
    score_file = os.path.join(output_dir, 'score.sc')
    cmd = [
        rosettapath,
        '-s {0}'.format(pdb),
        '-parser:protocol {0}'.format(xml),
        '-nstruct {0}'.format(nstruct),
        '-out:prefix {0}'.format(output_dir),
        '-out:file:score_only {0}'.format(score_file)
    ]
    subprocess.call(cmd)

    # Read in the results of scoring and return the mean
    # value over all runs
    time.sleep(5)
    df = pandas.read_csv(score_file, skiprows=1, sep='\s+')
    mean_score = df[filter_name].mean()

    return mean_score


def compute_protein_volume_metrics(pdb, vol_dir, path_to_ProteinVolume):
    """
    Use `ProteinVolume` to compute metrics on protein volume

    citation: Chen CR, Makhatadze GI (2015) ProteinVolume: calculating
    molecular van der Waals and void volumes in proteins. BMC
    Bioinformatics 16:1–6.

    Args:
        *pdb*: the path to a PDB file to analyze
        *vol_dir*: the path to a directory in which to store output
            from the program
        *path_to_ProteinVolume*: the path `ProteinVolume` script
    Returns:
        A tupple with the following elements in the following order:
            total_vol: the total volume of the protein
            void_vol: the total volume of voids
            vdw_vol: the total van der Waals volume
            packing_density: the packing density
    """

    # Make new directory for computing metrics on protein volume
    if not os.path.isdir(vol_dir):
        os.makedirs(vol_dir)

    # Copy PDB file to new directory
    new_pdb = os.path.join(
        vol_dir,
        os.path.basename(pdb)
    )
    subprocess.check_call(['cp', pdb, new_pdb])

    # Run command on new PDB file
    cmd = ' '.join([
        'java',
        '-jar',
        path_to_ProteinVolume,
        '"{0}"'.format(new_pdb)
    ])
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = process.communicate()

    # Read in results from the output
    out = out.decode('ascii')
    found_second_to_last_line = False
    for line in out.split('\n'):
        elements = line.split()
        if len(elements) == 0:
            continue
        if elements[0] == 'Protein':
            found_second_to_last_line = True
            continue
        if found_second_to_last_line:
            total_vol = elements[1]
            void_vol = elements[2]
            vdw_vol = elements[3]
            packing_density = elements[4]

    return(total_vol, void_vol, vdw_vol, packing_density)


def compute_abego_counts_in_loops(abego_string, dssp_string):
    """
    Compute counts for ABEGO 1-, 2-, and 3-mers in loops
    
    Args:
        *abego_string*: a string of the per-residue ABEGO types
            for each residue in the protein (upper case)
        *dssp_string*: a string of per-residue secondary structure
            for each residue in the protein (upper case; H=helix,
            E=strand, L=loop)
            
    Returns:
        A dictionary with all possible ABEGO 1-, 2-, and 3-mer
            sequences as keys and counts of these sequences in
            loops as values
            
    Code for doctest:
    >>> abego_string = 'AGA'
    >>> dssp_string = 'ELL'
    >>> d = compute_abego_counts_in_loops(abego_string, dssp_string)
    >>> (d['G'], d['A'], d['GA'])
    (1, 1, 1)
    >>> sum(d.values())
    3
    >>> abego_string = 'BAGBBEBBB'
    >>> dssp_string =  'LELLLHLLH'
    >>> d = compute_abego_counts_in_loops(abego_string, dssp_string)
    >>> (d['A'], d['B'], d['G'], d['GB'], d['BB'], d['GBB'])
    (0, 5, 1, 1, 2, 1)
    >>> sum(d.values())
    10
    """

    # Initiate a dictionary with all possible 1mer, 2mer, and 3mer ABEGO
    # sequences
    abego_1mers = list('ABEGO')
    abego_2mers = [
        ''.join(abego_2mer) for abego_2mer in
        list(itertools.product(abego_1mers, abego_1mers, repeat=1))
    ]
    abego_3mers = [
        ''.join(abego_3mer) for abego_3mer in
        list(itertools.product(abego_1mers, abego_2mers, repeat=1))
    ]
    abego_counts = {
        ''.join(list(abego_type)) : 0
        for abego_type in abego_1mers + abego_2mers + abego_3mers
    }

    # Go through sequence and record counts of 1mer, 2mer, and 3mer
    # ABEGOs
    abego_string = abego_string.upper()
    dssp_string = dssp_string.upper()
    assert len(abego_string) == len(dssp_string)
    for i in range(len(abego_string)):

        # Record data for 1mers if this position is in a loop
        if dssp_string[i] == 'L':
            abego_type_1mer = abego_string[i]
            abego_counts[abego_type_1mer] += 1
        else:
            continue

        # Record data for 2mers if the next position is in a loop
        if (i+1 >= len(abego_string)):
            continue
        if dssp_string[i+1] == 'L':
            abego_type_2mer = abego_string[i:i+2]
            abego_counts[abego_type_2mer] += 1
        else:
            continue

        # Record data for 3mers if the position after next is in a loop
        if (i+2 >= len(abego_string)):
            continue
        if dssp_string[i+2] == 'L':
            abego_type_3mer = abego_string[i:i+3]
            abego_counts[abego_type_3mer] += 1
    return abego_counts


def compute_total_charge_of_seq_subset(sequence, list_of_sites_ns):
    """
    Compute the total charge of a subset of an amino-acid sequence
    
    Args:
        *sequence*: amino-acid sequence (string)
        *list_of_site_ns*: list of site numbers (integers) that defines
            the subset of the sequence to analyze, indexed starting at 1
    Returns:
        The total charge of the sites listed (float)
        
    Code for doctest:
    >>> sequence = 'HAERKKD'
    >>> list_of_sites_ns = [1]
    >>> compute_total_charge_of_seq_subset(sequence, list_of_sites_ns)
    0.5
    >>> list_of_sites_ns = [1, 2, 3, 4, 5, 6, 7]
    >>> compute_total_charge_of_seq_subset(sequence, list_of_sites_ns)
    1.5
    """
    
    # Define amino-acid charges
    amino_acid_charges = {
        'E' : -1,
        'D' : -1,
        'R' : 1,
        'K' : 1,
        'H' : 0.5
    }
    
    # Compute the total charge of the indicated subset of sites
    # and return the result
    sequence = sequence.upper()
    total_charge = 0
    for site in list_of_sites_ns:
        aa = sequence[site-1]
        if aa in amino_acid_charges.keys():
            total_charge += amino_acid_charges[aa]
    return total_charge


def compute_per_residue_energies(pdb):
    """
    Compute per-residue energies for each term in the score function

    Args:
        *pdb*: the path to an input PDB file

    Returns:
        A dataframe with columns giving energies and rows giving
            residues
    """

    # Make a list of all score terms in the energy function
    score_terms = [
        str(score_term).replace('ScoreType.', '')
        for score_term in list(sf.get_nonzero_weighted_scoretypes())
    ]

    # Initiate a dictionary to store per-residue scores
    scores_dict = {
        key : []
        for key in ['res_n', 'res_aa', 'energy'] + score_terms
    }

    # Read in and score pose
    pose = pyrosetta.pose_from_pdb(pdb)
    sf(pose)

    # Make a dataframe with per-residue scores for each energy term
    # and the total score
    for res_n in list(range(1, pose.size()+1)):
        scores_dict['res_n'].append(res_n)
        scores_dict['res_aa'].append(pose.residue(res_n).name1())
        scores_dict['energy'].append(
            pose.energies().residue_total_energy(res_n)
        )
        for score_term in score_terms:
            scores_dict[score_term].append(
                pose.energies().residue_total_energies(res_n)[
                    pyrosetta.rosetta.core.scoring.score_type_from_name(
                        score_term
                    )
                ]
            )
    scores_df = pandas.DataFrame(scores_dict)
    return scores_df


if __name__ == '__main__':
    doctest.testmod()
