#!/bin/bash

cd /Users/calebbroodo/sos_analysis/

# Newer systematics:
#

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p3.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p4.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p5.txt input_table.txt
./getBusy.sh
#cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p5_TEST.txt input_table.txt
#./getBusy.sh
#
cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p6.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/DCA1p7.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/NhitsFit16.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/NhitsFit24.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/Vz12.txt input_table.txt
./getBusy.sh

cp -r /Users/calebbroodo/sos_analysis/tables_for_input/200GeV/Run11_v6_systematics/Vz21.txt input_table.txt
./getBusy.sh


