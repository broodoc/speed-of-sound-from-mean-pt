#!/bin/bash
echo ----====----====----====----====----====----====----====----====----====----

root -l -x <<EOF
.L Knee_Fitting.C
Knee_Fitting()
.q
EOF
echo ----====----====----====----====----====----====----====----====----====----

root -l knee_fitting.root -e 'new TBrowser'
