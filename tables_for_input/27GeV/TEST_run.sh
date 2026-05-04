#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

# Newer systematics:
#

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/27GeV/DCA2p0_TEST.txt input_table.txt
#cp -r /Users/calebbroodo/sos_analysis/tables_for_input/27GeV/DCA2p0_NoGlauber_TEST.txt input_table.txt
./getBusy.sh
