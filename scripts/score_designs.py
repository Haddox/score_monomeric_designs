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

# Import custom scripts
scriptsdir = os.path.dirname(__file__)
sys.path.append(scriptsdir)
import external_scoring

# Import and initialize `PyRosetta` and `Rosetta`
import pyrosetta
import pyrosetta.rosetta
import pyrosetta.distributed.io as io

# Score the PDB files
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
    
    # Conver the dictionary to a dataframe with the PDB file name as an index
    rs_scores = pandas.DataFrame.from_dict({pdb_file_name : rs_scores_dict}, orient='index')        
            
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

    # Remove the temporary results directory
    subprocess.check_call(['rm', '-r', resultsdir])

    # Merge scores
    assert all(rs_scores.index.values == external_scores.index.values)
    rs_metrics = list(rs_scores.columns.values)
    external_metrics = list(external_scores.columns.values)
    unique_external_metrics = [m for m in external_metrics if m not in rs_metrics]
    merged_scores = pandas.concat(
        [rs_scores, external_scores[unique_external_metrics]],
        axis=1
    )

    return merged_scores

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
pyrosetta.init(extra_options='-beta')
scorefxn = pyrosetta.get_score_function()

# Import the content of the input XML files, and store them in a dictionary
# of the form {metric name : script contents}
script_files = glob.glob(os.path.join('scripts/xmls', "*.xml"))
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


#---------------------------------------------------------------
# Use `Jug` to compute the scores for each PDB file in parallel
# and the results to a file
#---------------------------------------------------------------
# Score the PDB files
write_results(
    out_file_name,
    [ score_pdb_file(f, output_dir) for f in input_pdb_files ]
)