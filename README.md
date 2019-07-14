# Pipeline for computing an array of biophysical metrics for monomeric proteins

## Summary

This pipeline computes an array of biophysical metrics for an input set of monomeric proteins. The metrics include the ones included in the study [Rocklin et al., 2017, Science](http://science.sciencemag.org/content/357/6347/168), and much of the code was derived from that paper. There are also additional metrics that we added since then.

Note: There are many ways that this pipeline can be improved! Please feel free to make improvements that you come up with and push them to the repository. That could be anything from making better documentation, expanding the number of metrics compute, or fixing errors that cause the pipeline to crash.

## Organization of code

* `scripts/`: is a directory with a series of `Python` scripts that encode the scoring pipeline.
  * `sidechain_TdS_values.csv`: a file with estimates for losses in side-chain entropy upon folding, taken from Table 1 of Doig et al., 1995, Protein Science.
* `test_scoring.ipynb`: a notebook that runs the pipeline on ten structures from Rocklin et al., and then checks to see that the results from the pipeline match the results of the publication.
* `environment.yml`: a file listing a number of dependecies that are installable via `Conda` (see below).

## Installing external dependencies

Carrying out the pipeline requires multiple external dependencies. Unfortunately, the full set of required external dependencies are only available on the Baker lab server at this time. This includes:
* `PyRosetta`, as installed using `Conda` (see below).
* a few of the custom `Python` scripts for carying out the pipeline (e.g., `scripts/make_fragments.py` and `scripts/np_aa_burial.py`) refer to other scripts and programs on the Baker lab server that have not yet been extracted.

### Dependencies that are installable using `Conda`

Nearly all dependencies are encoded in the file called `environment.yml`. If you're working on the Baker lab server, these dependencies can all be installed using [`Conda`](https://conda.io/docs/index.html). To do so, first clone this repository. Then, in the root directory of the repository, execute the command:

    conda env create -f environment.yml -n {env_name}

where `env_name` is whatever name you'd like to call the environment.

### File with paths to external dependences that cannot be installed with `Conda`

Dependencies in `np_aa_burial.py`:

/software/rosetta/versions/v2019.01-dev60566/bin/rosetta_scripts.hdf5.linuxgccrelease

Dependencies in `make_fragments.py`:

fragment_tools = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/tools/fragment_tools/'
psipred = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/psipred3.21/'
scripts_dir = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/bakerlab_scripts/boinc/'
nnmake = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/nnmake/pNNMAKE.gnu'
csbuild = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/csbuild/'
cm_scripts = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/cm_scripts/bin/'
rosetta = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/main/source/bin/'

core_clusters, buried_np.xml



## How to run the pipeline

Running the pipeline involves two command-line arguments:

First, activate the `Conda` environment described above using the command:

    source activate {environment_name}

where `environment_name` is the name of the environment. This will give the pipeline access to many of the required external dependencies.

Second, use `jug` to execute the pipeline:

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

You can use `sbatch` to parallelize the jobs across multiple CPUs. For an example of a file that can be used as input for `sbatch`, see the file called `results/test_scoring/run.sbatch`, which is generated by `test_scoring.ipynb`. To submit the job to `slurm`, the notebook simply executes the command: `sbatch results/test_scoring/run.sbatch`.

## Checking the status of the run
You can also use `jug` to check the "status" of the job, i.e., how many jobs have and have not been completed:

    jug status --jugdir {jugdir} scripts/score_designs.py {path_to_folder_with_pdbs} {path_to_output_folder}

Unfortunately, this takes a while since it still needs to initialize all the Rosetta XML filters, which takes several minutes. I haven't figured out a way around this.

## What to do if the pipeline fails

Sometimes the pipeline fails in ways that are obvious, where the script raises a clear error (e.g., KeyError). In these cases, the code needs to be changed in order for it to work correctly. Of note: any temporary files and directories should be deleted before rerunning the pipeline (see below).

However, there are other times when the pipeline fails for reasons that are still mysterious to me. In these cases, I can delete temporary files and directories, rerun the pipeline immediately, and it works.

In each of the above cases, I suggest deleting all temporary files and directories before rerunning the pipeline. Below is a list of possible directories/files to check for:
  * there can be a variety of temporary files in the directory in which you executed the `jug` command from step 2 from above.
  * there can be `temp_*` directories in the output directory specified by `{path_to_output_folder}`.
  * `jug` will sometimes create `.lock` files in the directory `{jugdir}/locks/`. These should be deleted.

## Description of metrics

Most metrics are described in the supplemental material of [Rocklin, 2017, Science](http://science.sciencemag.org/content/357/6347/168). Below are descriptions of some of the new metrics:

**Metrics related to fragment quality:**
* `avg_all_frags_site_{N}`: for site `N` in the protein, the average RMS between the designed 9mer fragment centered in primary sequence upon that site and a set of 9mer fragments sampled from the PDB that are also centered in primary sequence upon that site.
* `avg_all_frags_in_H`, `avg_all_frags_in_E`, `avg_all_frags_in_L`: each of these metrics is a float that gives the average fragment quality of all 9mers centered on sites in a given secondary structure (H: helix, E: strand, L: loop). Specifically, each value is computed by considering all sites within a given secondary structure, and averaging the corresponding site-specific values given in the set of columns: `avg_all_frags_site_{N}`.

**Metrics related to per-residue Rosetta energies of fragments in primary sequence:**
* `energies_per_residue_site_{N}`: the total Rosetta energy for site `N` in the protein.
* `avg_per_residue_energies_{N}mer_{i}`: for a given fragment of size `N` starting at site `i` in the protein, the average per-residue total Rosetta energies of all sites in the fragment. Currently, this metric is computed for Ns of: 2, 3, 4, or 5.
* `avg_energy_for_{N}mers`, `min_energy_for_{N}mers`, `max_energy_for_{N}mers`: the average, minimum, and maximum value of all Nmers of a given size `N`. These values are computed from `avg_per_residue_energies_{N}mer_{i}`, across all sites `i`.

**Metrics related to per-residue energies of 3D neighborhoods:**
* `neighborhood_site_{i}_{D}A`: for site `i`, a comma-delimited string with integers giving the number of each residue in the protein that has a side-chain/side-chain contact with site `i`. Two residues are defined as contacting if there is at least one pair of atoms between residues that are within a distance of `D` Angstroms.
* `n_neighbors_site_{i}_{D}A`: for site `i`, the number of contacting residues in `neighborhood_site_{i}_{D}A`, including residue `i` itself.
* `avg_per_res_energy_of_site_{i}_neighborhood_{D}A`, `min_per_res_energy_of_site_{i}_neighborhood_{D}A`, `max_per_res_energy_of_site_{i}_neighborhood_{D}A`: for site `i`, the average, min, or max per-residue total Rosetta energy of all neighboring residues in 3D space, using a distance cutoff of `D` Angstroms for defining neighbors as described above for `neighborhood_site_{i}_{D}A`.
* `*_energy_of_{D}A_neighborhoods`: the average, min, or max of all values of `avg_per_res_energy_of_site_{i}_neighborhood_{D}A` across all sites `i` in the protein, using a distance cutoff of `D` Angstroms.
* `*_charge_of_{D}A_neighborhoods`: the average, min, or max charge of all neighborhoods across all sites `i` in the protein, using a distance cutoff of `D` Angstroms.


**Metrics related to counts of pairwise amino-acid contacts in 3D space**:
* `n_{i}{j}_3d_contacts_{D}A`: the number of times amino-acid `i` contacts amino-acid `j` in the structure. Residues are defined as contacting if they are neighbors in 3D space using a distance cutoff of `D` Angstroms for defining neighbors as described above for `neighborhood_site_{i}_{D}A`.

**Metrics related to buried hydrophobic surface area**:
* `buried_over_exposed_np_AFILMVWY`: `buried_np_AFILMVWY` / `exposed_np_AFILMVWY`

**Metrics related to side-chain entropy**:
* `{i}_TdS_{j}`: the change in side-chain entropy using TdS values from Table 1 of Doig et al., 1995, Protein Science, where `i` refers to the set of values from Table 1 (there are multiple estimates from different groups) and `j` is the layer in the protein (core, boundary, or surface) as defined using the side-chain-neighbor algorithm.

**Metrics related to packing and cavities**
* pack: measures packing density using the PackStat filter
* cavity_volume: uses the PackStat filter to estimate the total volume of voids in a structure
* holes: void volume computed using Holes filter
* `ProteinVolume_total_vol`: the total volume of the protein, including both the void volume and the van der Waals volume (see below) as computed using the `ProteinVolume` program (Chen CR, Makhatadze GI (2015) ProteinVolume: calculating molecular van der Waals and void volumes in proteins. BMC Bioinformatics 16:1–6.)
* `ProteinVolume_void_vol`: the void volume of the protein, computed using the `ProteinVolume` program
* `ProteinVolume_vdw_vol`: the van der Waals volume of the protein, computed using the `ProteinVolume` program.
* `ProteinVolume_packing_density`: the packing density of the protein (van der Waals volume / total volume), computed using the `ProteinVolume` program

**Metrics related to ABEGO types**
* `abego_counts_in_loops_*`: counts of 1-, 2-, and 3-mer ABEGO strings in loops for all possible strings
























***
