# trb_dirich_threshold
Program to set DiRICH thresholds at the HADES RICH

To compile the HADESthreshscan you need a version of ROOT (tested with >5.34), trbnet and boost-libraries.
The ROOT, boost and trbnet libraries should be filled in the Makefile: TRBNETDIR, ROOTDIR, BOOSTDIR.
The current setting of all of these varibles forsees a compilation on lxhadeb06 as hadaq.

Before using the program you need to set the correct LD_LIBRARY_PATH via ". setLD"

The usage of the program can be seen by using the --help /-h command
Examples of standard tasks are:

./HADESthreshscan_v1.C -b 0 -t 0 50
./HADESthreshscan_v1.C -l 0 -t 0 50
./HADESthreshscan_v1.C -f path/to/threshold.thr -l 0 -t 0 50
The -b/--baseline-scan performs a standard baselinescan. 
The option-parameter specifies which DiRICHes shall be scanned.
Here 0 specifies, that all DiRICHes shall be scanned.

The -t/--set-threshold sets the thresholds to a certain mV threshold above threshold. 
Here the first option-parameter (here 0) specifies the thresholds of which DiRICHes should be set (here 0 equals all DiRICHes). 
The second parameter equals the threshold to set in mV above threshold.

The -l/--load-baseline loads the baselines from a file which can be specified by the -f/--loading-file. 
If no file is specified the latest produced *.thr-file in the same directory is used. 
The option again defines which DiRICHes shall load baselines. 
Here 0 lets all DiRIChes load baselines if they exist in the specified file.
