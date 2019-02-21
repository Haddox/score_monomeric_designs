# Pipeline for computing an array of biophysical metrics for monomeric proteins

This pipeline computes an array of biophysical metrics for an input set of monomeric proteins. The metrics include the ones included in the study [Rocklin, 2017, Science](http://science.sciencemag.org/content/357/6347/168), and much of the code was derived from that paper. There are also additional metrics that we added since then.

There are many ways that this pipeline can be improved! Please feel free to make improvements that you come up with and push them to the repository. That could be anything from making better documentation, expanding the number of metrics compute, or fixing errors that cause the pipeline to crash.

## External dependencies

* Rosetta
* PyRosetta
* Scripts for fragment picking
* Everything in the `conda` environment

## Organization of code

* `scripts/`: is a directory with all code in the scoring pipeline
* `test_scoring.ipynb`: a notebook that runs the pipeline on ten structures from Rocklin et al., and then checks to see that the results from the pipeline match the results of the publication.

## How to run the pipeline

Running the pipeline involves two command-line arguments:

1) Activate the `conda` environment described above using the command:

    source activate {environment_name}

where `environment_name` is the name of the environment. This will give the pipeline access to many of the required external dependencies.

2) Use `jug` to execute the pipeline:

    jug execute --jugdir {jugdir} scripts/score_designs.py {path_to_directory_with_input_pdbs} {path_to_output_directory}

* Inputs for the above command:
  * `{jugdir}`: the path to a directory where `jug` will store a variety of files needed to execute and keep track of the run.
  * `{path_to_directory_with_input_pdbs}`: the path to a directory with all input PDBs. Currently, all input PDBs need to be in this single directory.
  * `{path_to_output_directory}`: the path to a directory where all output will be stored, including temporary directories and files used in the computation and the final scores file.

* Outputs of the above command:
  * `scores.csv`: a CSV file with all computational metrics for each design, where columns are different metrics and rows are different designs. This file will be written to the directory specified by the input variable called `{path_to_output_directory}`.
  * the pipeline also generates a large number of temporary directories and files for each input structure, which are automatically deleted if the pipeline runs to completion. However, if the pipeline fails in the middle of the run, some of these files may not get deleted. See the below section called `What to do if the pipeline fails for some of the designs`

## Example Jupyter notebook that runs the pipeline with a handful of test PDBs

The notebook called `test_scoring.ipynb` runs the pipeline on ten structures from Rocklin et al., and then checks to see that the results from the pipeline match the results of the publication.

## How to run the pipeline using `sbatch`

You can use `sbatch` to parallelize the jobs across multiple CPUs.

You can also use `jug` to check the "status" of the job, i.e., how many jobs have and have not been completed:

    jug execute --jugdir {jugdir} scripts/score_designs.py {path_to_folder_with_pdbs} {path_to_output_folder}

Unfortunately, this takes a while since it still needs to initialize all the Rosetta XML filters, which takes several minutes. I haven't figured out a way around this.

## What to do if the pipeline fails

Sometimes the pipeline fails in ways that are obvious, where the script raises a clear error (e.g., KeyError). In these cases, the code needs to be changed in order for it to work correctly. Of note: any temporary files and directories should be deleted before rerunning the pipeline (see below).

However, there are other times when the pipeline fails for reasons that are still mysterious to me. In these cases, I can delete temporary files and directories, rerun the pipeline immediately, and it works.

In each of the above cases, I suggest deleting all temporary files and directories before rerunning the pipeline. Below is a list of possible directories/files to check for:
  * there can be a variety of temporary files in the directory in which you executed the `jug` command from step 2 from above.
  * there can be `temp_*` directories in the output directory specified by `{path_to_output_folder}`.
  * `jug` will sometimes create `.lock` files in the directory `{jugdir}/locks/`. These should be deleted.


## Description of metrics

Most metrics are described in the supplemental material of [Rocklin, 2017, Science](http://science.sciencemag.org/content/357/6347/168).

* `avg_all_frags_per_site`: a string with comma-delimited entries that report the average fragment quality within a sliding window centered upon each site. Specifically, these values are reported for each site that sits at the center of a 9mer window in the protein. Values are ordered by sites, and there are no values for the first four or last four sites in the protein since these sites do not sit at the center of any 9mer windows. Fragment quality is quantified as the RMS between the design and 9mer fragment. The reported number of each site is the average RMS over all fragments in a given window.
* `avg_all_frags_in_H`, `avg_all_frags_in_E`, `avg_all_frags_in_L`: each of these metrics is a float that gives the average fragment quality of all 9mers centered on sites in a given secondary structure (H: helix, E: strand, L: loop). Specifically, each value is computed using by considering all sites within a given structure, and averaging the corresponding site-specific values given in `avg_all_frags_per_site`.

## Ways to improve pipeline

* figure out why pipeline will sometimes fail with large batches of designs, but then works if it's rerun


## Missing filters:
* `cavity_volume`
    * Currently gives rise to the error:
        Error: Element 'CavityVolume': This element is not expected. Expected is one of ( AASynthesisFitnessCost, AlaScan, AlignmentAAFinder, AlignmentGapInserter, AngleToVector, AtomCount, AtomicContact, AtomicContactCount, AtomicDistance, AverageDegree ).
    * Parisa said that this may not yet be part of master
