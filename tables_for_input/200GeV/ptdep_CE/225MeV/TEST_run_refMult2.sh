#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

# Newer systematics:
#

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/ptdep_CE/225MeV/DCA1p8_refMult2.txt input_table.txt
./getBusy.sh
