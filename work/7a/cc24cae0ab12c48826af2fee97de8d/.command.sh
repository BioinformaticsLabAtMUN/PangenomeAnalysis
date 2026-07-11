#!/bin/bash -ue
cd-hit -i consolidated.faa \
    -o cdhit_output.faa \
    -c 0.65 \
    -n 4 \
    -aL 0.75 \
    -G 1 \
    -g 1 \
    -T 12 \
    -M 10000
