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
    external_scores.set_index('Unnamed: 0', inplace=True)

    # Remove the temporary results directory and the dssp file
    subprocess.check_call(['rm', '-r', resultsdir])
    subprocess.check_call(['rm', '-f', pdb_file_name+'.dssp'])

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

    # Sum the `rama_prepro` and `p_aa_pp` terms and make columns that indicate whether
    # structures have helices and/or strands.
    scores_df['rama_prepro_and_p_aa_pp'] = scores_df['rama_prepro'] + scores_df['p_aa_pp']
    scores_df['has_helices'] = scores_df['secondary_structure'].apply(lambda x: 'H' in x)
    scores_df['has_strands'] = scores_df['secondary_structure'].apply(lambda x: 'E' in x)

    # Compute the frequency of each amino acid across the entire protein, in
    # alpha helices, in beta strands, or in different layers.
    amino_acids = IUPAC.IUPACProtein.letters
    amino_acids = list(amino_acids)
    assert len(amino_acids) == 20
    for aa in amino_acids:

        # compute amino-acid frequencies across the entire protein
        scores_df['freq_{0}'.format(aa)] = scores_df['sequence'].apply(lambda x: x.count(aa)/len(x))

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
        lambda row: 100.0 * sum([row['freq_{0}'.format(aa)] for aa in hydrophobic_amino_acids]),
        axis=1
    )
    for layer in layers:
        scores_df['percent_hydrophobic_AFILMVWY_{0}'.format(layer)] = scores_df.apply(
            lambda row: 100.0 * sum([row['{0}_freq_{1}'.format(layer, aa)] for aa in hydrophobic_amino_acids]),
            axis=1
        )

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

# Import the content of the input XML files, and store them in a dictionary
# of the form {metric name : script contents}
script_files = glob.glob(os.path.join(scriptsdir, 'xmls', "*.xml"))
script_files = script_files[:1]
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
input_pdb_files = input_pdb_files = glob.glob(os.path.join(pdb_dir, '*.pdb'))
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
