#!/bin/bash

qsub -q short.q -o ./log-bash -e ./log-bash -N tests -b y "matlab -nodesktop -r 'addpath ./test; test_basic_functions_publish;'"

