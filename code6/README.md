# Antenna Radiation Pattern Simulation

## Overview
This project contains a C++ program that simulates the radiation pattern of a phased array antenna using different windowing techniques. The program calculates the antenna gain as a function of the angle \( \theta \) and generates a graphical output using Gnuplot.

## Features
- Computes the far-field radiation pattern for a phased array antenna.
- Supports different windowing techniques for element weighting:
  - **Blackman**
  - **Chebyshev**
  - **Taylor**
- Generates a Gnuplot script to visualize the radiation pattern in dB.
- Allows customization of the number of elements, spacing, frequency, and phase shift.

## Dependencies
Ensure you have the following installed before running the program:
- A C++ compiler (e.g., `g++`)
- Gnuplot (for visualization)

## Compilation
To compile the program, use the following command:

```bash
 g++ -o antenna_simulation antenna_simulation.cpp -std=c++11 -O2 -lm
```

## Usage
Run the executable with the following parameters:

```bash
 ./antenna_simulation <N> <L> <frequency> <delta_phi> <method> <ripple_or_sll>
```

### Parameters:
- `N` : Number of elements in the array.
- `L` : Total length of the antenna array.
- `frequency` : Operating frequency (Hz).
- `delta_phi` : Phase shift between elements (radians).
- `method` : Windowing method (`blackman`, `chebyshev`, or `taylor`).
- `ripple_or_sll` : Ripple level for Chebyshev or sidelobe level for Taylor (dB).

### Example:
```bash
 ./antenna_simulation 40 4 1.2e9 0.1 taylor -30.0
```

## Output
- A data file `data_cartesian.dat` containing the computed radiation pattern.
- A Gnuplot script `cartesian_plot.gnu` to generate the plot.
- A PNG image `cartesian_plot.png` displaying the radiation pattern.

## Notes
- The program assumes a far-field condition for the calculations.
- If Gnuplot is not installed, the script will generate the data file but won't plot the results automatically.

## License
This project is released under the MIT License. Feel free to modify and use it for your research or projects.

