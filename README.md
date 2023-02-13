# NLFLeigenvectors

Instructions to run tests:

* go to the tests directory - you will find 
  * a SAVE for h-BN 2D
  * example-LR directory - this is a test at the frequency of the first exciton
  * example-SHG directory - this is a test at the frequency of the first peak in Chi^(2)
* inside each example directory you will find
  * run.sh - submission script for Lobster to be used with yambo-devel-nl executables (edit path)
  * copy this and the SAVE directory to Lobster, and submit run.sh - it will:
    * execute command.sh, which creates new input files with slightly different total simulation times
    * calculate the Coulomb integrals
    * run yambo_nl and ypp_nl in a loop

* In order to execute the analysis you will need to
  * retrieve the job_0# directories from Lobster
  * have the latter and SAVE/ns.db1 in the directory where you intend to run the python processing
  * run process_all_kpts.py

