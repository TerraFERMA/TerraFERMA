#!/bin/bash
###########
# quick and dirty bash script to set the environment and run the test harness with N threads
###########
Nthreads=1
. ./environment
tfsimulationharness --test -n $Nthreads compare_solvers.shml
