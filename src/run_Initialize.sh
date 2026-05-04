#!/bin/bash

# Compile

root -l <<EOF
.L StParameters.cxx+
.L Initiate_input.cpp+
Initiate_input();
.q
EOF
