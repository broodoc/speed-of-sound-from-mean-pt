#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p8_refMult2.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p8.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/54GeV/DCA2p0_TEST.txt input_table.txt
./getBusy.sh



