#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

# Newer systematics:
#

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p5_TEST.txt input_table.txt
./getBusy.sh
