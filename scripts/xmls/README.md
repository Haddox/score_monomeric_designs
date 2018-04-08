This directory contains a set of RosettaScripts used to score designs using the same filters as [Rocklin et al., 2017, Science](http://science.sciencemag.org/content/357/6347/168). There is one script for each filter for ease of testing. Specifically, files are named `*.xml`, where the star indicates the `name="*"` argument for the filter, as specified in the script. We make sure that these scripts are working as expected in [`test_score_protocols.ipynb`](../../test_score_protocols.ipynb).


Here is a list of all the filters from the Rocklin et al. study. I am currently in the process of making a single XML for each. I will track my progress by putting a `#` at the beginning of the line of each filter for which I've made its own XML.

#  <Add filter_name="res_count_core_SCN" />
#  <Add filter_name="res_count_core_SASA" />
#  <Add filter_name="percent_core_SCN" />
#  <Add filter_name="percent_core_SASA" />
#    <Add filter_name="contact_all" />
#    <Add filter_name="contact_core_SCN" />
#    <Add filter_name="contact_core_SASA" />
#  <Add filter_name="degree" />
#  <Add filter_name="entropy" />
  <Add filter_name="dslf_quality_check" /> `will leave out for now if our designs don't have disulfides`
  <Add filter_name="mean_dslf" /> `will leave out for now if our designs don't have disulfides`


  <Add filter_name="cavity_volume" /> #`having difficulty`
#  <Add filter_name="ss_sc" />
#  <Add filter_name="helix_sc" />
#  <Add filter_name="loop_sc" />
  <Add filter_name="exposed_total" /> `RuntimeError: EXCN_Base::what()`
  <Add filter_name="exposed_hydrophobics" /> `RuntimeError: EXCN_Base::what()`
  <Add filter_name="exposed_polars" /> `RuntimeError: EXCN_Base::what()`
  <Add filter_name="fxn_exposed_is_np" /> `RuntimeError: EXCN_Base::what()`
  <Add filter_name="holes" /> `ERROR: Value of inactive option accessed: -holes:dalphaball`; Gabe provided this path of a file using a command-line arg
#  <Add filter_name="bb" />

  <Add filter_name="buried_np" /> `same problem as exposed filters`
  <Add filter_name="buried_over_exposed" /> `same problem as exposed filters`
  <Add filter_name="buried_minus_exposed" /> `same problem as exposed filters`
#  <Add filter_name="pack" />
  <Add filter_name="mismatch_probability" /> `working, but with non-local files`
    Add filter="rotamer_boltz_core_avg" />
    Add filter_name="degree_core_SCN" />
    Add filter_name="degree_core_SASA" />
#    <Add filter_name="unsat_hbond" />
    <Add filter_name="unsat_hbond2" /> `Doesn't recognize filter?`
#    <Add filter_name="one_core_each" />
#    <Add filter_name="two_core_each" />
#    <Add filter_name="ss_contributes_core" />
#    <Add filter_name="AlaCount" />
    Add filter_name="faulty_fragments_tolerant" />
    Add filter_name="faulty_fragments" />

Need nres?
