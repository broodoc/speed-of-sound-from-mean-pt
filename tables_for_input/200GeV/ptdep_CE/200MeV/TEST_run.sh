#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

# Newer systematics:
#

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/ptdep_CE/200MeV/DCA2p0.txt input_table.txt
./getBusy.sh
