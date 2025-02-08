# Antenna Radiation Pattern Simulation

## Overview
This program simulates the radiation pattern of an antenna array by computing the electric field at various angles. The results are saved in a data file and plotted using Gnuplot.

## Features
- Models an antenna array with multiple emitters
- Computes received power based on emitter distance
- Generates a radiation pattern in dB
- Produces a visual representation using Gnuplot

## Requirements
- C++ compiler (GCC, Clang, MSVC, etc.)
- Gnuplot (for visualization)

## Compilation
To compile the program, run:
```sh
 g++ -o antenna_simulation antenna_simulation.cpp -O2 -std=c++11
```

## Usage
Run the compiled executable:
```sh
 ./antenna_simulation
```
The program will generate:
- `data_cartesian.dat`: Contains computed radiation pattern data
- `cartesian_plot.png`: The resulting radiation pattern plot

## Expected Output
Upon execution, the program:
1. Computes power at different angles
2. Saves the results in `data_cartesian.dat`
3. Generates a Gnuplot script `cartesian_plot.gnu`
4. Runs Gnuplot to create `cartesian_plot.png`

The output plot shows the gain in dB versus the angle in degrees, illustrating the radiation pattern of the antenna array.

## Notes
- Ensure Gnuplot is installed and available in your system path.
- Adjust parameters like the number of emitters and spacing within the source code to explore different configurations.
