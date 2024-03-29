# Code and results
This directory contains code and results from "A covariance-enhanced approach to
multi-tissue joint eQTL mapping with application to
transcriptome-wide association studies"

MTeQTLResults contains both the data and example code for extracting estimated eQTL weights. These weights can be used, for example, to perform MultiXcan-type testing. 

MutliXcanResults contains the p-values for each gene tested in Section 5 of the main manuscript. 

Simulations contains an R script which was used for all simulation studies in the main manuscript. To recreate these simulation results, the Simulations.sh bash script should be executed -- this will run 7500 separate jobs: one for each replicate in the simulation studies. The output will be saved in a directory (which must be created) named Results. With the results in hand, the plots from Fig 2 can be recreated exactly. Of course, one will have to be careful to modify filepaths correctly according to their own computing environment.
