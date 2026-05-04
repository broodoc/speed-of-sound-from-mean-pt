THESE .root files are used for statistical bootstrapping
for error bar estimation (on charged differential yield).
Statistical bootstrap is done by splitting the full dataset
into 10 equally distributed sets of .root files (and merging
those .root files).
Note: Run in the following order - 
1. ls output*root | sort -R > all_files.txt
2. create_ten_groups.sh
3. merge_ten.sh
4. Statistical bootstrap code
Note: The RCF code and configuration is the same as Dec_8_2024.
