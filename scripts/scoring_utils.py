"""
Modules for scoring protein structures using `PyRosetta`
"""

import pandas
import doctest
import pyrosetta
import pyrosetta.rosetta
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from pyrosetta.rosetta.core.select.residue_selector import LayerSelector
from pyrosetta.bindings.utility import bind_method

@bind_method(HBondSet)
def __iter__(self):
    for i in range(1, self.nhbonds() + 1):
        yield self.hbond(i)

@bind_method(HBondSet)
def hbonds(self):
    return [hb for hb in self]

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

if __name__ == '__main__':
    doctest.testmod()
