# Antenna Simulation - Fortran

## Overview
This Fortran program simulates the radiation pattern of an antenna array and computes the received power at various angles. It models an antenna array and calculates the field intensity based on parameters such as the number of elements, frequency, spacing, and phase shift.

The program generates a data file containing the computed gain values in dB and creates a Gnuplot script to visualize the results.

## Features
- Computes the average power received at a given range and angle.
- Models an antenna array with user-defined parameters.
- Normalizes and converts power values to dB.
- Generates Gnuplot-compatible data for visualization.

## Requirements
- A Fortran compiler (e.g., `gfortran`)
- Gnuplot for visualization

## Compilation & Execution
### Compilation
To compile the program, use:
```sh
 gfortran -o antenna_simulation antenna_simulation.f90
```

### Execution
Run the program with the following arguments:
```sh
 ./antenna_simulation <num_elements> <frequency> <d> <radius> <delta_phi>
```
Where:
- `<num_elements>`: Number of elements in the antenna array
- `<frequency>`: Operating frequency (Hz)
- `<d>`: Element spacing (meters)
- `<radius>`: Distance at which power is computed (meters)
- `<delta_phi>`: Phase shift between elements (radians)

### Example
```sh
 ./antenna_simulation 40 1.2e9 0.5 2000 0.1
```

## Output
- The computed power data is stored in `data_cartesian2.dat`.
- A Gnuplot script `cartesian_plot2.gnu` is generated.
- The program automatically runs Gnuplot to produce `cartesian_plot2.png`.

## Visualization
To manually generate the plot, run:
```sh
 gnuplot cartesian_plot2.gnu
```
This will produce a PNG file showing the antenna radiation pattern.

## Notes
- Ensure Gnuplot is installed to visualize the results.
- The program computes power in ideal conditions without interference or external noise.
