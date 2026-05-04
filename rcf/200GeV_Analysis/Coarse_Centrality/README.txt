The set of code written in this directory /star/u/cbroodo/200GeV_Analysis/RCF/200GeV_Analysis/Coarse_Centrality
is purpose-built for analyzing mean pT fluctuations using the traditional RefMultCorr binning for centrality
defined in Run16 by https://drupal.star.bnl.gov/STAR/system/files/Run16_centralityDef.pdf. The relevant
centrality regime is the 0-5%, 5-10%, in this context. The analysis in this directory is used to support 
the argument that finer binning is required to analyze a measureable slope parameter in the ultra-central 
regime.

11.13.24 - Global StRefMultCorr does not work. Ask Jan Vanek to send you Github link to download
his StRefMultCorr geared for Run14 + Run16 analysis and work with it accordingly. Upload the 
StRefMultCorr/ to StRoot/StSoundVelocity and find a way to link it/reference it in your 
analysis (ask Jan if you need help). 

11.27.24 - StRefMultCorr works now (included Jan's local version) but difficulty arises
when trying to implement WeightingHisto. Disk quota exception.

11.28.24 - Checking to make sure code procures files w/o implementing weighting histo.
If they work, will store files in local /Desktop/200GeV_Run16/RCF_Files/Nov_28_2024/

4.11.25 - Trying to match the Nch distributions btw Run16 and Run11

Optimization notes:
Would be useful to create an enum class for the different configurations of centrality vs analysis 
determinations 
