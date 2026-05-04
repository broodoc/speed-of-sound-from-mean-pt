#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

# Newer systematics:
#

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/54GeV/ptdep/NoCut/DCA2p0_TEST.txt input_table.txt
./getBusy.sh
