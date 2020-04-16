#!/bin/bash

gnuplot -p -e "set format y '%.1e'; set log y; plot '<&3' w lp" 3< <(grep deltaT log.slmMeltPoolFoam | awk '{print $3}')
