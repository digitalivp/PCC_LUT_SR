#!/bin/bash


## Folders
TMC13DIR=../mpeg-pcc-tmc13/
RESULTS=test/experiment/
COND=octree-raht/lossy-geom-lossy-attrs/

## Input parameters
FILE=longdress_vox10_1300;
RATE=6;

OUT_FOLDER="${TMC13DIR}${RESULTS}${COND}${FILE}/r0${RATE}/";

## LUT based SR
echo ----- Running SR ------
MFILE=mfile.m
echo "addpath('./scripts/');" > ${MFILE}
echo "sr_lut_frac_core('"$FILE"', "$RATE", '"$OUT_FOLDER"');" >> ${MFILE}
octave-cli ${MFILE}
rm ${MFILE}
