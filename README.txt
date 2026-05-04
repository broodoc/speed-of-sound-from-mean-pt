This directory streamlines the speed of sound analysis
(should require minimal human input, though there's
no gaurantee that errors won't occur). The human input
that is required:

SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP 
Copy + paste this in terminal: 

file /usr/local/Cellar/root/6.32.08/lib/root/libCore.so
arch -x86_64 zsh
uname -m

SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP 

1. Alter input_table.txt (look at all elements to ensure
you are making the correct references)

2. run: ./getBusy.sh
This runs through the entire analysis in the following steps:
 - creates prepped_files/name_from_input/ directory (this stores all of the 
resultant output, including the final graphs_output.root file)
 - isolates QA plots
 - centrality determination (to change centrality specs based 
on centrality estimator, alter Parameterization.h file)
 - statistical bootstrapping
 - analysis (spectra fitting, mean pT and dN/deta calculation, 
sos fitting, generating final graphs) 
 - saves input_table.txt and terminal output log and all 
relevant files to prepped_files/name_from_input/

3. Navigate to the relevant directory and examine the recently created
.root file (root -l .root -> new TBrowser) to see SOS analysis

Optional: To rerun selected parts of the analysis,
refer to the individual *sh files that correspond
to the steps above in the src/ directory 

Optional: To free up space from multiple output root
And log files, simply run the ./cleanup.sh script
That gets copied into your prepped files directory
When running ./getBusy.sh

graph_output.root file structure:

fOut
 ├── pT study
 │    ├── Absolute quantities
 │    ├── Normalized (to ref class) quantities
 │    │    ├── pT-Integrated Nch
 │    │    │    ├── fit_range_0.0-0.5_GeV
 │    │    │    ├── fit_range_0.0-0.6_GeV
 │    │    │    └── ...
 │    └── Referenced Nch
 ├── spectra study
 └── general


---------------------------------------------------------------------------

1. Input table (.txt file with directory reference w/ RCF files, refMult name, pTrefMult name, prepped_files directory name, Config information - eta range for analysis, eta and pT range for centrality estimator, Run, collision system (COL or FXT), collision energy) -> create prepped_files directory name/ -> store input_table.txt in prepped_files directory name/ -> DONE
2.  Check that directory reference w/ RCF files exists -> Store directory reference to RCF files in prepped_files directory name/ -> DONE
3.  Remove 0 byte files -> DONE
4.  Check that refMult name, pTrefMult name exists by checking the first output*.root -> DONE
5.  Generate “cleaned_output.*” isolating relevant histos (based on refMult name, pTrefMult name) -> DONE
6.  Hadd output files to temp -> DONE
7.  (Optional) Isolate QA plots -> convert to beautiful pdfs and store in prepped_files directory name/
8.  Centrality determination of temp -> accept relevant directory path name and file name, refMult name, as input argument for FindPercentiles.cpp and store in prepped_files directory name/ -> DONE
9.  Copy merging code into directory reference w/ RCF files and run ./run_merging -> DONE
10.  Statistical bootstrapping -> fix to relevant directory path name (as input argument for Bootstrap.C and run statistical bootstrapping) and run ./run_Bootstrapping.sh -> store yield_errors in prepped_files directory name/ -> DONE
11.  Copy yield_errors into Analysis_Eta_pT_comp -> DONE
12.  Analysis -> relevant directory path name, refMult name, pT_refMult name as input arguments for Run_Analysis -> store graph_output.root in prepped_files directory name/ -> Store fitlog in prepped_files directory name/ -> DONE
13.  Save terminal output to a analysis_log_(date and time).txt -> Store in prepped_files directory name/ -> DONE
14.  Merge relevant files (Trajectum and CMS) -> DONE
15.  Time the run time (estimate the amount of runtime) -> DONE
