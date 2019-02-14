# Pipeline for computing an array of biophysical metrics for monomeric proteins

## External dependencies

* PyRosetta
* Rosetta
* Scripts for fragment picking
* Everything in the conda environment

## How to run the pipeline



## Organization of code

* `scripts/` : is a directory containing all code in the scoring pipeline
* `data/` : contains data for analysis

## Explanation of metrics

* `avg_all_frags_per_site`: a string with comma-delimited entries that report the average fragment quality across a sliding window. Specifically, these values are reported for each site that sits at the center of a 9mer window in the protein. Values are ordered by sites, and there are no values for the first four or last four sites in the protein since these sites do not sit at the center of any 9mer windows. Fragment quality is quantified as the RMS between the design and 9mer fragment. The reported number of each site is the average RMS over all fragments in a given window.
* `avg_all_frags_in_H`, `avg_all_frags_in_E`, `avg_all_frags_in_L`: each of these metrics is a float that gives the average fragment quality of all 9mers centered on sites in a given secondary structure (H: helix, E: strand, L: loop). Specifically, each value is computed using by considering all sites within a given structure, and averaging the corresponding site-specific values given in `avg_all_frags_per_site`.

## Missing filters:
* `cavity_volume`
    * Currently gives rise to the error:
        Error: Element 'CavityVolume': This element is not expected. Expected is one of ( AASynthesisFitnessCost, AlaScan, AlignmentAAFinder, AlignmentGapInserter, AngleToVector, AtomCount, AtomicContact, AtomicContactCount, AtomicDistance, AverageDegree ).
    * Parisa said that this may not yet be part of master



## To do:
* Ask Will about the `cavity_volume` filter
