"""
Python script for computing a variety of structural metrics that, for the most
part are NOT encoded in RosettaScripts XMLs
"""

import unittest
import pandas
import subprocess
import glob
import os
from os import path
import sys
import json

import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.score as score

#TODO refactor to add informative top-level module name for `scripts` directory
script_dir = path.dirname(__file__)
sys.path.append(script_dir)
#from test_data.rocklin_rd4_examples import test_pdb_files
assert sys.path.pop() == script_dir


def generate_enhanced_score_summary_table(
        input_poses,
        scriptsdir,
        resultsdir,
        output_score_file_prefix
    ):
    """
    Generate a table with RosettaScripts and non-RosettaScripts scores.

    For each input pose, this script first uses Rosetta to compute scores,
    including those specified in the input XML files. Next, this script
    computes additional scores using the following custom Python scripts,
    which are external dependencies:
        `./scripts/make_fragments.py`
        `./scripts/enhance_scorefile.py`

    Args:
        input_poses - { string pdb path : pyrosetta.distributed.packed_pose.PackedPose }

        scriptsdir - path of a directory with `make_fragments.py` and
            `enhance_scorefile.py`.

        resultsdir - path where results will be stored

        output_score_file_prefix - the prefix to the output score files
            in the path of `resultsdir`

    Returns:
        A CSV with scores computed from ONLY from the input XML files for each
            pose, which as the path: `{resultsdir}/{output_score_file_prefix}.csv`.
            This file doesn't have the additional scoring metrics computed using
            the custom Python scripts.

        For each pose, a directory with information on the fragments
            genreated for each pose, with the path `{resultsdir}/
            fragment_qual_results/*_fragments/`, where `*` is the
            name of a given pose

        A file that summarizes all fragment quality scores with the path:
            {resultsdir}/score.sc.frag_quals

        An enhanced score file generated by `enhance_scorefile.py`, which has
            the path `{resultsdir}/{output_score_file_prefix}_enhanced.csv`.
    """

    if not path.exists(resultsdir):
        os.makedirs(resultsdir)

    # Compute scores for each pose using the input XML RosettaScripts
    print("Generating input scorefile", file=sys.stderr)

    scoring = score.ScorePoseTask()

    input_table = pandas.DataFrame.from_dict(
        { name : packed_pose.to_dict(scoring(p)) for name, p in input_poses.items() },
        orient="index"
    )

    del input_table["pickled_pose"]

    # Write the results to an output file as input for `enhance_scorefile.py`,
    # adding an additional column that gives the complete path to each pdb file,
    # minus the PDB extension, so that it is compatible with `enhance_scorefile.py`
    input_table['description'] = input_table.index.map(lambda x: os.path.splitext(x)[0])
    csv_inputfile = os.path.join(resultsdir, '{0}.csv'.format(output_score_file_prefix))
    input_table.to_csv(csv_inputfile, index='False')

    fragment_resultsdir = path.join(resultsdir, "fragment_qual_results")
    if not path.exists(fragment_resultsdir):
        os.makedirs(fragment_resultsdir)

    # Generate and score fragments with the script `./scripts/make_fragments.py`
    print("Generating and scoring fragments with `make_fragments.py`", file=sys.stderr)
    print("Storing the results in the directory: {0}".format(fragment_resultsdir))
    input_pdb_list = list(input_poses.keys())
    if scriptsdir[-1] == '/':
        scriptsdir = scriptsdir[:-1]
    subprocess.check_call(
        ['python', '{0}/make_fragments.py'.format(path.relpath(scriptsdir, fragment_resultsdir)), '-pdbs'] + list(map(path.abspath, input_pdb_list)),
        stderr=subprocess.PIPE,
        cwd=fragment_resultsdir,
    )

    # Concatenate the `frag_qual.summary` files for each input pose into a single
    # file as input for `enhance_scorefile.py`
    print("Concatenating fragment quality scores", file=sys.stderr)
    frag_qual_summary_files = glob.glob(
            path.join(fragment_resultsdir, '*_fragments/frag_qual.summary'))
    master_frag_qual_file = '{0}.frag_quals'.format(csv_inputfile)
    with open(master_frag_qual_file, 'w') as f:
        subprocess.check_call(
            ['cat'] + frag_qual_summary_files,
            stderr=subprocess.PIPE, stdout=f
        )

    # Use `./scripts/enhance_scores.py` to compute additional non-RosettaScripts
    # scores
    print("Computing additional non-RosettaScripts scores with `enhance_scorefile.py` with the input CSV file {0}".format(csv_inputfile), file=sys.stderr)
    result_data = subprocess.check_output(
        ['python2', '{0}/enhance_scorefile.py'.format(scriptsdir), csv_inputfile])

    print("result_data", result_data)
    return json.loads(result_data.decode())

class TestExternalScoring(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.input_poses = {
            f : io.pose_from_file(f)
            for f in test_pdb_files[:2]
        }

    def test_external_scoring(self):
        result = generate_enhanced_score_summary_table(
            self.input_poses,
            path.dirname(__file__),
            path.join(path.dirname(__file__), "../results/test_external_scoring_results"),
            "test"
        )

        self.assertTrue(len(result) == len(self.input_poses))

        for r in result:
            self.assertTrue("avg_all_frags" in r)

if __name__ == "__main__":
    unittest.main()
