#!/bin/bash

file /usr/local/Cellar/root/6.32.08/lib/root/libCore.so

arch -x86_64 zsh <<EOF

uname -m

source /usr/local/Cellar/root/6.32.08/bin/thisroot.sh

make clean

make

EOF
