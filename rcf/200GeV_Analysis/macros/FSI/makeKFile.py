import os

os.system("root -b -q 'FSI.C++(\"3.0\")'")
os.system("root -b -q 'FSI.C++(\"3.5\")'")
os.system("root -b -q 'FSI.C++(\"4.0\")'")
os.system("root -b -q 'FSI.C++(\"4.5\")'")
os.system("root -b -q 'FSI.C++(\"5.0\")'")
os.system("root -b -q 'FSI.C++(\"5.5\")'")
os.system("root -b -q 'FSI.C++(\"6.0\")'")
os.system("root -b -q 'FSI.C++(\"6.5\")'")
os.system("root -b -q 'FSI.C++(\"7.0\")'")
os.system("root -b -q 'FSI.C++(\"7.5\")'")
os.system("root -b -q 'FSI.C++(\"8.0\")'")
os.system("root -b -q 'FSI.C++(\"8.5\")'")
os.system("root -b -q 'FSI.C++(\"9.0\")'")
os.system("root -b -q 'FSI.C++(\"9.5\")'")
os.system("root -b -q 'FSI.C++(\"10.0\")'")

# R VALUES NEED TO BE MANUALLY in mergeKFiles.C !!!!!!
os.system("root -b -q 'mergeKFiles.C({\"output_fsi_3.0fm.root\",\"output_fsi_3.5fm.root\",\"output_fsi_4.0fm.root\",\"output_fsi_4.5fm.root\",\"output_fsi_5.0fm.root\",\"output_fsi_5.5fm.root\",\"output_fsi_6.0fm.root\",\"output_fsi_6.5fm.root\",\"output_fsi_7.0fm.root\",\"output_fsi_7.5fm.root\",\"output_fsi_8.0fm.root\",\"output_fsi_8.5fm.root\",\"output_fsi_9.0fm.root\",\"output_fsi_9.5fm.root\",\"output_fsi_10.0fm.root\"})'")

