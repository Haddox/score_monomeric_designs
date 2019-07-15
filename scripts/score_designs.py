"""
Score designs using `PyRosetta` and a series of input RosettaScript XMLs
"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import pandas
import re
import scipy.stats
import numpy
import math
import jug
import subprocess
from Bio.Alphabet import IUPAC
import itertools

# Import custom scripts
scriptsdir = os.path.dirname(__file__)
sys.path.append(scriptsdir)
import external_scoring
import scoring_utils

# Import and initialize `PyRosetta` and `Rosetta`
import pyrosetta
import pyrosetta.rosetta
import pyrosetta.distributed.io as io


# Score the PDB files
@jug.TaskGenerator
def score_pdb_file(pdb_file_name, output_dir):
    """
    Score an input PDB file using an array of Rosetta-based metrics

    PyRosetta computes these scores using the metrics encoded in the
    RosettaScripts XML files in the directory `./scripts/xmls/`.

    Args:
        `pdb_file_name`: The path to a PDB file (string)
        `output_dir`: The path to a directory where temporary files
            will be stored and then deleted.

    Returns:
        A dataframe with scores for the PDB file
    """

    #------------------------------------------------------------------------
    # Read in the PDB as a pose and use PyRosetta to compute sequence- and
    # structure-based metrics, including many used in Rocklin et al., 2017,
    # Science
    #------------------------------------------------------------------------

    # Read in the PDB file and convert it to a pose
    pose = pyrosetta.pose_from_pdb(pdb_file_name)

    # Compute scores from energies from the score function and store
    # the scores in a dictionary with the form {score name : value}
    rs_scores_dict = {}
    scorefxn(pose)
    for score in scorefxn.get_nonzero_weighted_scoretypes():
        re_pattern = re.compile(r'ScoreType\.(?P<score_name>\w+)')
        score_name = re_pattern.search(str(score)).group('score_name')
        rs_scores_dict[score_name] = pose.energies().total_energies()[score]

    # Compute XML-based scores and add these scores to the dictionary
    # with all the scores
    for (filter_name, pr_filter) in pr_filters.items():

        # See if the scoring works
        try:
            print("Scoring the design {0} with the filter: {1}".format(pdb_file_name, filter_name))
            rs_scores_dict[filter_name] = pr_filter.score(pose)

        # If there is a RuntimeError, put "Error" for the given pose/feature
        except RuntimeError:
            #logger.exception("Error processing feature: %s" % feature_name)
            rs_scores_dict[filter_name] = 'RuntimeError'

    # Add the pose's sequence and secondary structure to the score dictionary
    sequence = pose.sequence()
    rs_scores_dict['sequence'] = sequence
    DSSP.apply(pose)
    rs_scores_dict['secondary_structure'] = pose.secstruct()

    # Add lists of bools specifying whether residues are in the core, surface,
    # or boundary based on the side-chain neighbor algorithm, storing the bools
    # as a string of zeros and ones
    core_bools = core_selector.apply(pose)
    boundary_bools = boundary_selector.apply(pose)
    surface_bools = surface_selector.apply(pose)
    rs_scores_dict['core_bools'] = ''.join(
        numpy.array(core_bools).astype(int).astype(str)
    )
    rs_scores_dict['boundary_bools'] = ''.join(
        numpy.array(boundary_bools).astype(int).astype(str)
    )
    rs_scores_dict['surface_bools'] = ''.join(
        numpy.array(surface_bools).astype(int).astype(str)
    )

    # Add values quantifying the number of residues in each layer
    rs_scores_dict['nres_core'] = sum(core_bools)
    rs_scores_dict['nres_boundary'] = sum(boundary_bools)
    rs_scores_dict['nres_surface'] = sum(surface_bools)

    # Make sure that each residue is assigned to one and only one layer
    for (i,j,k) in zip(
        rs_scores_dict['core_bools'],
        rs_scores_dict['boundary_bools'],
        rs_scores_dict['surface_bools']
        ):
        assert(sum([int(n) for n in [i,j,k]]) == 1)

    # Conver the dictionary to a dataframe with the PDB file name as an index
    rs_scores = pandas.DataFrame.from_dict({pdb_file_name : rs_scores_dict}, orient='index')

    # Compute additional scores based on the Rosetta output
    if 'buried_np' in rs_scores.columns.values:
        rs_scores['buried_np_per_res'] = rs_scores['buried_np'] / rs_scores['nres']
        rs_scores['buried_over_exposed'] = rs_scores['buried_np'] / (rs_scores['exposed_hydrophobics'] + 0.01)
        rs_scores['buried_minus_exposed'] = rs_scores['buried_np'] - rs_scores['exposed_hydrophobics']


    #------------------------------------------------------------------------
    # Use `Python` scripts from Rocklin et al., 2017, Science to score each
    # input PDB using additional sequence- and structure-based metrics
    #------------------------------------------------------------------------
    # Compute external scores
    resultsdir  = os.path.join(
        output_dir,
        'temp_{0}/'.format(
            os.path.splitext(pdb_file_name)[0].replace('/', '_')
        )
    )
    output_score_file_prefix = 'score'
    external_scores_dict = external_scoring.generate_enhanced_score_summary_table(
        {pdb_file_name : io.pose_from_file(pdb_file_name)},
        scriptsdir,
        resultsdir,
        output_score_file_prefix
    )
    external_scores = pandas.DataFrame.from_dict(external_scores_dict)
    external_scores.rename(columns={'Unnamed: 0':'pdb_path'}, inplace=True)
    external_scores.set_index('pdb_path', inplace=True)

    # Merge scores
    assert all(rs_scores.index.values == external_scores.index.values)
    rs_metrics = list(rs_scores.columns.values)
    external_metrics = list(external_scores.columns.values)
    unique_external_metrics = [m for m in external_metrics if m not in rs_metrics]
    scores_df = pandas.concat(
        [rs_scores, external_scores[unique_external_metrics]],
        axis=1
    )


    #------------------------------------------------------------------------
    # Compute additional metrics related to amino-acid frequencies in
    # different parts of the structure
    #------------------------------------------------------------------------

    # Sum the `rama_prepro` and `p_aa_pp` terms and make columns that indicate
    # whether structures have helices and/or strands.
    scores_df['rama_prepro_and_p_aa_pp'] = \
        scores_df['rama_prepro'] + scores_df['p_aa_pp']
    scores_df['has_helices'] = \
        scores_df['secondary_structure'].apply(lambda x: 'H' in x)
    scores_df['has_strands'] = \
        scores_df['secondary_structure'].apply(lambda x: 'E' in x)

    # Compute the frequency of each amino acid across the entire protein, in
    # alpha helices, in beta strands, or in different layers.
    amino_acids = IUPAC.IUPACProtein.letters
    amino_acids = list(amino_acids)
    assert len(amino_acids) == 20
    for aa in amino_acids:

        # compute amino-acid frequencies across the entire protein
        scores_df['freq_{0}'.format(aa)] = \
            scores_df['sequence'].apply(lambda x: x.count(aa)/len(x))

        # ... only in alpha helices or beta strands
        for (ss, ss_id) in [('helices', 'H'), ('strands', 'E')]:
            scores_df['{0}_freq_{1}'.format(ss, aa)] = scores_df.apply(
                lambda row: scoring_utils.ComputeAminoAcidFrequencyInResidueSubset(
                    aa_seq = row['sequence'],
                    per_residue_ids = row['secondary_structure'],
                    aa = aa,
                    subset_id = ss_id
                ),
                axis=1
            )

        # ... only in the core, boundary, or surface
        layers = ['core', 'boundary', 'surface']
        for layer in layers:
            scores_df['{0}_freq_{1}'.format(layer, aa)] = scores_df.apply(
                lambda row: scoring_utils.ComputeAminoAcidFrequencyInResidueSubset(
                    aa_seq = row['sequence'],
                    per_residue_ids = row['{0}_bools'.format(layer)],
                    aa = aa,
                    subset_id = '1'
                ),
                axis=1
            )

    # Compute the frequency of each dipeptide across the entire protein
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            dipeptide = amino_acids[i] + amino_acids[j]
            scores_df['freq_{0}'.format(dipeptide)] = \
                scores_df['sequence'].apply(
                    lambda x: x.count(dipeptide) / float(len(x) - 1)
                )

    # Compute joint frequency of hydrophobic amino acids defined by the list below
    hydrophobic_amino_acids = list('AFILMVWY')
    scores_df['percent_hydrophobic_AFILMVWY'] = scores_df.apply(
        lambda row: 100.0 * sum(
            [row['freq_{0}'.format(aa)] for aa in hydrophobic_amino_acids]
        ),
        axis=1
    )
    for layer in layers:
        scores_df['percent_hydrophobic_AFILMVWY_{0}'.format(layer)] = \
            scores_df.apply(
                lambda row: 100.0 * sum([
                    row['{0}_freq_{1}'.format(layer, aa)]
                    for aa in hydrophobic_amino_acids
                ]), axis=1
        )

    #------------------------------------------------------------------------
    # Compute additional metrics related to buried hydrogen bonds
    #------------------------------------------------------------------------

    # Compute the number of buried hydrogen bonds using VSASA to determine burial,
    # with a burial cutoff of 0.1, as is used in the `BuriedUnsatHbonds` filter
    vsasa_burial_cutoff = 0.1
    buried_h_bonds = scoring_utils.find_buried_hbonds(pose)
    if buried_h_bonds:
        scores_df['buried_hbonds_info'] = ';'.join(buried_h_bonds)
        scores_df['n_buried_hbonds'] = len(buried_h_bonds)
    else:
        scores_df['buried_hbonds_info'] = buried_h_bonds
        scores_df['n_buried_hbonds'] = 0

    #------------------------------------------------------------------------
    # Compute additional metrics related to per-site fragment quality
    #------------------------------------------------------------------------

    # Compute the average fragment quality among all 9mers
    # centered upon a given site. To so do, I will first
    # read in data on site-specific 9mer fragment quality
    fragment_resultsdir = os.path.join(
        resultsdir,
        "fragment_qual_results"
    )
    frag_qual_file = glob.glob(os.path.join(
        fragment_resultsdir,
        '*_fragments',
        'frag_qual.dat'
    ))
    assert len(frag_qual_file) == 1, frag_qual_file
    frag_qual_file = frag_qual_file[0]
    names = ['fragment_size', 'position', 'fragment_index', 'rms', 'x', 'y']
    frag_qual_df = pandas.read_csv(
        frag_qual_file, sep='\s+', names=names
    )
    assert sum(frag_qual_df['fragment_size'] == 9) == len(frag_qual_df), \
        "Not all fragments are 9mers"

    # Next, for each site in the protein, I will compute the average
    # fragment quality (i.e., RMS between design and fragment) of all
    # 9mer fragments that are centered on that site. Thus, there will
    # be no data for the first four and last four residues, since
    # they are never in the middle of any fragments.
    frag_qual_df['centered_position'] = frag_qual_df['position'] + 4
    avg_frag_qual_dict = {
        key : []
        for key in ['centered_position', 'avg_rms']
    }
    for position in list(set(frag_qual_df['centered_position'])):
        avg_frag_qual_dict['centered_position'].append(position)
        avg_rms = frag_qual_df[
            frag_qual_df['centered_position'] == position
        ]['rms'].mean()
        avg_frag_qual_dict['avg_rms'].append(avg_rms)
        scores_df['avg_all_frags_site_{0}'.format(position)] = avg_rms
    avg_frag_qual_df = pandas.DataFrame.from_dict(avg_frag_qual_dict)
    avg_frag_qual_df.sort_values(
        'centered_position', ascending=True, inplace=True
    )
    avg_all_frags_per_site = list(avg_frag_qual_df['avg_rms'])

    # Next, compute the average fragment quality over all sites
    # in different elements of secondary structure, specifically:
    # helices, strands, and loops
    assert len(scores_df) == 1, scores_df
    secondary_structure = list(scores_df['secondary_structure'][0])[4:-4]
    ss_types = ['L', 'H', 'E']
    frag_qual_ss_dict = {
        key : []
        for key in ss_types
    }
    for (ss, avg_frag_qual) in zip(secondary_structure, avg_all_frags_per_site):
        frag_qual_ss_dict[ss].append(avg_frag_qual)
    for ss_type in ss_types:
        if len(frag_qual_ss_dict[ss_type]) == 0:
            scores_df['avg_all_frags_in_{0}'.format(ss_type)] = numpy.nan
        else:
            scores_df['avg_all_frags_in_{0}'.format(ss_type)] = \
                pandas.Series(frag_qual_ss_dict[ss_type]).mean()

    #------------------------------------------------------------------------
    # Compute additional metrics related to per-site Rosetta energy in local
    # neighborhoods in primary sequence and 3D neighborhoods, as well as other
    # metrics related to 3D structural contacts
    #------------------------------------------------------------------------

    # Make a dataframe with the total energy of each residue, and record
    # each residue's energy in the main dataframe returned at the end of this
    # function
    nres = pose.total_residue()
    energies_df = scoring_utils.compute_per_residue_energies(pdb_file_name)
    energies_df.sort_values('res_n', inplace=True, ascending=True)
    energies_df.set_index('res_n', inplace=True)
    for (i, row) in energies_df.iterrows():
        scores_df['energy_per_residue_site_{0}'.format(i)] = row['energy']

    # Make sure the per-residue energies add up to approximately
    # the total energy
    assert (sum(energies_df['energy']) - scorefxn(pose)) < 1e-3
        
    # For each score term, compute summary statistics on per-residue
    # energies, including the total energy, as well as individual score
    # terms.
    for energy_col in energies_df.columns.values:
        if energy_col in ['res_aa']:
            continue        
        
        # First, compute the min, max, mean, and standard deviation of
        # each term. Also compute the mean of the 5 worst or best scores.
        sorted_vals = energies_df[energy_col].sort_values(ascending=True)
        scores_df['min_per_site_{0}'.format(energy_col)] = sorted_vals.min()
        scores_df['max_per_site_{0}'.format(energy_col)] = sorted_vals.max()
        scores_df['mean_per_site_{0}'.format(energy_col)] = sorted_vals.mean()
        scores_df['std_per_site_{0}'.format(energy_col)] = sorted_vals.std()
        scores_df['avg_lowest_per_site_{0}'.format(energy_col)] = \
            sorted_vals[:5].mean()
        scores_df['avg_highest_per_site_{0}'.format(energy_col)] = \
            sorted_vals[-5:].mean()

        # Next, compute the energy of all n-mers in the structure and then
        # compute summary statistics from that list
        fragment_sizes = [2, 3, 4, 5]
        for fragment_size in fragment_sizes:
            fragment_energies = []
            for res_n in range(1, nres-fragment_size+2):
                fragment_n = [res_n + i for i in range(fragment_size)]
                avg_per_residue_energy_fragment_i = \
                    sum(energies_df.loc[fragment_n][energy_col]) / fragment_size
                fragment_energies.append(avg_per_residue_energy_fragment_i)
                scores_df['avg_per_residue_{0}_{1}mer_{2}'.format(
                    energy_col, fragment_size, res_n
                )] = avg_per_residue_energy_fragment_i
            scores_df['avg_{0}_for_{1}mers'.format(energy_col, fragment_size)] = \
                numpy.mean(fragment_energies)
            scores_df['min_{0}_for_{1}mers'.format(energy_col, fragment_size)] = \
                numpy.min(fragment_energies)
            scores_df['max_{0}_for_{1}mers'.format(energy_col, fragment_size)] = \
                numpy.max(fragment_energies)
            scores_df['std_{0}_for_{1}mers'.format(energy_col, fragment_size)] = \
                numpy.std(fragment_energies)

    # For each residue in the protein, get a list of all neighboring residues
    # where both residues have at least one side-chain atom below the specified
    # distance cutoff
    distance_cutoffs = [3, 5]
    for distance_cutoff in distance_cutoffs:

        # Make a list of neighbors
        energies_df['neighborhood_{0}'.format(distance_cutoff)] = \
            energies_df.apply(
                lambda row: scoring_utils.get_residue_neighbors(
                    pdb_file_name, row.name, distance_cutoff
                ),
                axis=1
            )

        # For each residue's neighborhood, compute the number of residues in the
        # neighborhood, including the residue used to define the neighborhood
        energies_df['n_neighbors_{0}'.format(distance_cutoff)] = \
            energies_df['neighborhood_{0}'.format(distance_cutoff)].apply(
                lambda x: len(x)
            )

        # ... then compute the average per-residue Rosetta total energy of all
        # residues in the neighborhood
        energies_df['energy_of_neighborhood_{0}'.format(distance_cutoff)] = \
            energies_df.apply(
                lambda row: sum(
                    energies_df.loc[
                        row['neighborhood_{0}'.format(distance_cutoff)]
                    ]['energy']) / row['n_neighbors_{0}'.format(distance_cutoff)],
                axis=1
            )

        # ... then compute the total charge of all residues in the
        # neighborhood
        energies_df['charge_of_neighborhood_{0}'.format(distance_cutoff)] = \
            energies_df.apply(
                lambda row: scoring_utils.compute_total_charge_of_seq_subset(
                    sequence,
                    row['neighborhood_{0}'.format(distance_cutoff)]
                ), axis=1
            )
        
        # ... then compute the average distance in primary sequence between the
        # residue used to define the neighborhood and each of its neighbors
        energies_df[
            'avg_dist_in_primary_sequence_to_neighbors_{0}'.format(distance_cutoff)
        ] = \
            energies_df.apply(lambda row: numpy.mean([
                abs(row.name - neighbor_n)
                for neighbor_n in row['neighborhood_{0}'.format(distance_cutoff)]
                if neighbor_n != row.name
            ]), axis=1)

    # Add the above data into the main dataframe returned at the
    # end of this function
    for distance_cutoff in distance_cutoffs:

        # First, add site-specific data
        for (i, row) in energies_df.iterrows():
            scores_df[
                'neighborhood_site_{0}_{1}A'.format(row.name, distance_cutoff)
            ] = \
                ','.join(map(str, row['neighborhood_{0}'.format(distance_cutoff)]))
            scores_df[
                'avg_per_res_energy_of_site_{0}_neighborhood_{1}A'.format(
                    row.name, distance_cutoff
                )
            ] = \
                row['energy_of_neighborhood_{0}'.format(distance_cutoff)]
            scores_df[
                'n_neighbors_site_{0}_{1}A'.format(row.name, distance_cutoff)
            ] = \
                row['n_neighbors_{0}'.format(distance_cutoff)]
            scores_df[
                'avg_dist_in_primary_sequence_to_neighbors_site_{0}_{1}A'.format(
                    row.name, distance_cutoff
                )
            ] = \
                row[
                'avg_dist_in_primary_sequence_to_neighbors_{0}'.format(
                    distance_cutoff
                    )
                ]

        # Then, add summary statistics that take all sites into account
        scores_df['avg_energy_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['energy_of_neighborhood_{0}'.format(distance_cutoff)].mean()
        scores_df['min_energy_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['energy_of_neighborhood_{0}'.format(distance_cutoff)].min()
        scores_df['max_energy_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['energy_of_neighborhood_{0}'.format(distance_cutoff)].max()
        scores_df['std_energy_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['energy_of_neighborhood_{0}'.format(distance_cutoff)].std()
        
        scores_df['mean_charge_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['charge_of_neighborhood_{0}'.format(distance_cutoff)].mean()
        scores_df['min_charge_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['charge_of_neighborhood_{0}'.format(distance_cutoff)].min()
        scores_df['max_charge_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['charge_of_neighborhood_{0}'.format(distance_cutoff)].max()
        scores_df['std_charge_of_{0}A_neighborhoods'.format(distance_cutoff)] = \
            energies_df['charge_of_neighborhood_{0}'.format(distance_cutoff)].std()

        
    # Compute the number of pairswise 3D contacts for all amino-acid pairs
    # for each distance cutoff
    for distance_cutoff in distance_cutoffs:

        # First, initiate a dictionary with keys giving each pairwise
        # combination of amino acids, with each key being specified by two
        # single-letter codes in alphabetical order, and with values
        # set equal to zero
        aa_pairwise_contacts_dict = {}
        for (i, aa_i) in enumerate(amino_acids):
            for (j, aa_j) in enumerate(amino_acids):
                aa_pairwise_contacts_dict[tuple(sorted([aa_i, aa_j]))] = 0

        # Then, for each target site, iterate over all neighboring sites and
        # add each contact to the above dictionary, making sure to ignore
        # the target site when it comes up in the list of neighboring sites
        for (i, row) in energies_df.iterrows():
            target_aa = pose.residue(row.name).name1()
            neighboring_aas = [
                pose.residue(res_n).name1() for res_n in
                row['neighborhood_{0}'.format(distance_cutoff)]
                if res_n != row.name
            ]
            for neighboring_aa in neighboring_aas:
                aa_pairwise_contacts_dict[
                    tuple(sorted([target_aa, neighboring_aa]))
                ] += 1

        # Add the above counts to the main dataframe that is returned at the
        # end of this function
        for (aa_i, aa_j) in aa_pairwise_contacts_dict.keys():
            scores_df['n_{0}{1}_3d_contacts_{2}A'.format(
                aa_i, aa_j, distance_cutoff
            )] = aa_pairwise_contacts_dict[(aa_i, aa_j)]


    #------------------------------------------------------------------------
    # Compute metrics related to hydrophobic and polar surface area
    #------------------------------------------------------------------------
    scores_df['buried_over_exposed_np_AFILMVWY'] = \
        scores_df['buried_np_AFILMVWY'] / \
            scores_df['exposed_np_AFILMVWY']
    scores_df['buried_minus_exposed_np_AFILMVWY'] = \
        scores_df['buried_np_AFILMVWY'] - \
            scores_df['exposed_np_AFILMVWY']


    #------------------------------------------------------------------------
    # Compute metrics related to conformational entropy
    #------------------------------------------------------------------------
    # Read in estimates for the loss in conformational entropy upon residue
    # burial
    tds_df = pandas.read_csv('scripts/sidechain_TdS_values.csv')
    tds_df.rename(columns={'single-letter code': 'aa'}, inplace=True)
    assert len(tds_df['aa']) == len(set(tds_df['aa']))
    tds_groups = [
        'Picket_Sternberg', 'Abagyan_Totrov', 'Koehl_Delarue', 'Lee'
    ]

    # Compute loss in entropy for each layer in the protein, defining layers using
    # the side-chain-neighbors algorithm
    for layer in ['core', 'boundary', 'surface']:
        for tds_group in tds_groups:

            # First, compute TdS over all amino acids in a layer by multiplying
            # the TdS and frequency of each amino acid and summing these values
            # over all amino acids
            scores_df['{0}_TdS_{1}_per_res'.format(tds_group, layer)] = \
                scores_df.apply(
                    lambda row: sum([
                        float(row['{0}_freq_{1}'.format(layer, aa)]) * \
                            float(tds_df[tds_df['aa'] == aa][tds_group])
                        for aa in amino_acids
                    ]), axis = 1
                )

            # The above values are computed with amino-acid frequencies. Below,
            # I multiply the above TdS value by the number of residues in a layer
            # to get a total TdS value, NOT normalized by the number of residues
            scores_df['{0}_TdS_{1}'.format(tds_group, layer)] = \
                scores_df['{0}_TdS_{1}_per_res'.format(tds_group, layer)] * \
                    scores_df['nres_{0}'.format(layer)]


    #------------------------------------------------------------------------
    # Compute metrics using RosettaScripts command-line arguments. I will do this
    # for filters that for some reason aren't working with PyRosetta
    #------------------------------------------------------------------------
    xmls_to_run_on_commandline = glob.glob(
        os.path.join(scriptsdir, 'xmls_to_run_on_commandline/*.xml')
    )
    nstruct = 10
    for xml in xmls_to_run_on_commandline:
        filter_name = os.path.basename(xml).replace('.xml', '')
        output_dir = os.path.join(resultsdir, filter_name) + '/'
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        scores_df[filter_name] = scoring_utils.score_design_from_command_line(
            rosettapath,
            pdb_file_name,
            xml,
            nstruct,
            output_dir,
            filter_name
        )

    #------------------------------------------------------------------------
    # Add metrics on ABEGO types
    #------------------------------------------------------------------------
    # First, add a column that gives the per-residue ABEGO type of all sites
    abego_f = os.path.join(resultsdir, 'score.csv.abego_terms')
    with open(abego_f) as f:
        lines = f.readlines()
        assert len(lines) == 1, abego_f
        entries = lines[0].strip().split()
        assert len(entries) == 5, abego_f
        scores_df['abego_types'] = entries[2].upper()

    # Then, add counts of ABEGO-type strings for all 1-mer, 2-mer, and 3-mers
    # within loops
    abego_string = scores_df.iloc[0]['abego_types'].upper()
    dssp_string = scores_df.iloc[0]['dssp'].upper()
    assert len(abego_string) == len(dssp_string) == len(sequence), \
        "{0} vs. {1} vs. {2}".format(
        len(abego_string), len(dssp_string), len(sequence)
    )
    abego_counts_dict = scoring_utils.compute_abego_counts_in_loops(
        abego_string, dssp_string
    )
    for (abego_string, counts) in abego_counts_dict.items():
        scores_df['abego_counts_in_loops_{0}'.format(abego_string)] = counts

    # Then, compute counts of each amino acid in each ABEGO type in loops.
    # To do so, first initiate a dictionary of all possible aa-abego combos...
    abego_types = list('ABEGO')
    abego_aa_counts_in_loops_dict = {
        '{0}_{1}'.format(aa, abego) : 0
        for (aa, abego) in list(itertools.product(amino_acids, abego_types))
    }

    # ... then count the number of observed combos in loops
    for (abego, aa, ss) in zip(abego_string, sequence, dssp_string):
        if ss.upper() == 'L':
            abego_aa_counts_in_loops_dict['{0}_{1}'.format(aa, abego)] += 1
    for (combo, counts) in abego_aa_counts_in_loops_dict.items():
        scores_df['{0}_counts_in_loops'.format(combo)] = counts
        
    #------------------------------------------------------------------------
    # Add metrics on protein volume
    #------------------------------------------------------------------------
    vol_dir = os.path.join(resultsdir, 'protein_vol/')
    path_to_ProteinVolume = 'scripts/ProteinVolume_1.3/ProteinVolume_1.3.jar'
    vol_output = scoring_utils.compute_protein_volume_metrics(
        pdb_file_name, vol_dir, path_to_ProteinVolume
    )
    assert len(vol_output) == 4, len(vol_output)
    scores_df['ProteinVolume_total_vol'] = vol_output[0]
    scores_df['ProteinVolume_void_vol'] = vol_output[1]
    scores_df['ProteinVolume_vdw_vol'] = vol_output[2]
    scores_df['ProteinVolume_packing_density'] = vol_output[3]

    #------------------------------------------------------------------------
    # Compute per-residue values for a subset of metrics
    #------------------------------------------------------------------------
    per_res_metrics = [
        'exposed_np_AFILMVWY', 'buried_minus_exposed_np_AFILMVWY',
        'buried_over_exposed_np_AFILMVWY'
    ]
    for metric in per_res_metrics:
        scores_df['{0}_per_res'.format(metric)] = \
            scores_df[metric] / scores_df['nres']


    #------------------------------------------------------------------------
    # Remove temporary files and folders and return all data in the `scores_df`
    # dataframe
    #------------------------------------------------------------------------
    # Remove the temporary results directory and the dssp file
    subprocess.check_call(['rm', '-r', resultsdir])
    subprocess.check_call(['rm', '-f', pdb_file_name+'.dssp'])

    # Return all data in dataframe
    return scores_df

@jug.TaskGenerator
def write_results(results_file, sub_results):
    # Merge the scores of all PDBs into a single dataframe
    scores_all_pdbs = pandas.concat(sub_results)

    # Write the results to an outfile
    if not os.path.isdir(os.path.dirname(results_file)):
        os.makedirs(os.path.dirname(results_file))
    print("Writing the results to the file: {0}".format(results_file))
    scores_all_pdbs.to_csv(results_file)

#---------------------------------------------------------------
# Initialize relevant Movers and Filters in PyRosetta
#---------------------------------------------------------------
# Initialize a specific score function
pyrosetta.init(extra_options='-beta -holes:dalphaball /work/brunette/scripts/DAlphaBall.gcc')
scorefxn = pyrosetta.get_score_function()

# Tweak the score function so that it includes energies from
# backbone hydrogen bonds when computing per-residue energies
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(scorefxn)

# Import the content of the input XML files, and store them in a dictionary
# of the form {metric name : script contents}
script_files = glob.glob(os.path.join(scriptsdir, 'xmls', "*.xml"))
def _read(f):
    with open(f) as inf:
        return inf.read()
analysis_scripts = {
    os.path.splitext(os.path.basename(x))[0] : _read(x)
    for x in script_files
}

# Initialize the filters by reading in and processing the XML file for each filter
pr_filters = {

    filter_name : pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(script_contents)
        .get_filter(filter_name)
    for (filter_name, script_contents) in list(analysis_scripts.items())

}
DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
core_selector = scoring_utils.layer_selector_mover('core')
boundary_selector = scoring_utils.layer_selector_mover('boundary')
surface_selector = scoring_utils.layer_selector_mover('surface')

# Define path to Rosetta. TODO: make this a command-line argument
rosettapath = '/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease'

#---------------------------------------------------------------
# Read in command-line arguments and input PDB files
#---------------------------------------------------------------
# Read in command-line arguments using `argparse`
parser = argparse.ArgumentParser()
parser.add_argument("pdb_dir", help="the path to a directory with input PDB files")
parser.add_argument("output_dir", help="the path to a directory where the output will be stored")

args = parser.parse_args()
(pdb_dir, output_dir) = (args.pdb_dir, args.output_dir)
out_file_name = os.path.join(output_dir, 'scores.csv')

# Get the input PDB files
input_pdb_files = glob.glob(os.path.join(pdb_dir, '*.pdb'))
if len(input_pdb_files) == 0:
    raise ValueError("Could not find any PDB files in the directory: {0}".format(pdb_dir))


#---------------------------------------------------------------
# Use `Jug` to compute the scores for each PDB file in parallel
# and the results to a file
#---------------------------------------------------------------
# Score the PDB files
write_results(
    out_file_name,
    [ score_pdb_file(f, output_dir) for f in input_pdb_files ]
)
