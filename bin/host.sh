#!/bin/sh
#
# This script guesses the system id string in the same way configure does
# Provided for use in scripts that have something to do with PSI3
#
# Center for Computational Quantum Chemistry,
# University of Georgia,
# Athens, GA
# February 2000

#
# First try config.guess
# config.guess is free software provided with GNU autoconf
#
script_dir=`dirname $0`
if [ -x $script_dir/config.guess ]; then
   target=`$script_dir/config.guess`
fi
if [ -x $script_dir/config.local ]; then
   target=`$script_dir/config.local $target`
fi

echo $target
