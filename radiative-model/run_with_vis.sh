#!/bin/sh

echo -n "script running\n"

# compile the model and call the executable radiative_surface_balance
gcc -o radiative_surface_balance radiative_surface_balance.c -lm

# run the model
./radiative_surface_balance

# plot the results (if pylab is installed)
python plot.py
