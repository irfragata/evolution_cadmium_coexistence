This folder contains all the code necessary to run the analyses

If you want to repeat the main analyses run the Main_analyses.R script. The other scripts are either called by this script, or are used to perform parameter exploration and methods comparison.

- Main_analyses.R - This file has the whole code to run the main analyses of the paper. If the parameter evaluation is TRUE the code starts by calling the Running_cxr_pooledData.R and Running_cxr_Replicates.R scripts to estimate the model parameters. Then it provides all the steps to perform the analyses described in the manuscript. At the end of the script there is a call to the Figures.R script to obtain all figures in the manuscript (except Fig S2). Note that the code produces several intermediate files that are needed for figures. These data files are available in the Analyses folder in case you want to run the Figures script without running this one. 
Html files are available in the html folder, allowing to see the code and the output.

Warining: These two scripts take a long time to run.
- Running_cxr_pooledData.R - code to estimate the parameters for the Riker model using pooled data. This code takes a while to run. Intermediate files with the data produced in this script are available in the Analyses folder in case you are only interested in the data analyses. 
- Running_cxr_Replicates.R - code to estimate the parameters for the Riker model using data from each individual replicate. This code takes a long time to run. Intermediate files with the data produced in this script are available in the Analyses folder in case you are only interested in the data analyses.

- Figures.R - This file has the whole code to obtain the figures produced for the paper. The code needs several intermediate files that are produced in the main analyses. These data files are available in the Analyses folder in case you want to run the Figures script without running the main analyes. Html files are available in the html folder, allowing to see the code and the output.

The following files are to perform additional analyses

- Function_riker.R - contains the functions needed to run the methods D and E used in the Method Comparison script. 

- Method_Comparison.R - WARNING: this code takes a while to run! Code to test which methods show better performance when estimating parameters from the data. This code is used to produce figure S2.

- Parameter_exploration.R - Warning: this code takes a while to run! Code to do an initial grid search to test for the initial conditions with which to start the cxr run, to avoid encountering local subobtima likelihood peaks. To cut a bit the time that it takes to run, you can load the data from the Parameter_exploration.RData file.


- Session_Info - provides all the information of the packages used in the analyses



