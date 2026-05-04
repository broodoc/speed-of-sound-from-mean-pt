#!/bin/bash

for i in $(seq 1 10); do
    hadd -f merged_${i}.root $(cat group_${i}.txt)
done
