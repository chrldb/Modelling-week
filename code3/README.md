# Radiation Pattern Simulation

## Description
This program simulates the radiation pattern of an antenna array and calculates the necessary phase shifts to scan a given angular range.

## Features
- Computes the far-field radiation pattern for different phase shifts.
- Uses Gnuplot to generate graphical representations of the results.
- Determines the number of phase shifts (aims) required for scanning a given range.
- Outputs the power distribution in dB over different angles.

## Requirements
- A C++ compiler (e.g., `g++`)
- Gnuplot installed for visualization

## Compilation
To compile the program, use the following command:
```bash
 g++ -o radiation_pattern radiation_pattern.cpp -std=c++11 -lm
```

## Usage
To run the program, use the following command:
```bash
./radiation_pattern <N> <L> <frequency>
```
where:
- `<N>` is the number of antenna elements.
- `<L>` is the total length of the antenna array.
- `<frequency>` is the operating frequency in Hz.

## Example
For an array with 40 elements, a length of 4 meters, and a frequency of 1.2 GHz:
```bash
./radiation_pattern 40 4 1.2e9
```

## Output
- The number of phase shifts (aims) required.
- A data file (`radiation_pattern.dat`) containing computed power values.
- A plot (`radiation_pattern.png`) displaying the radiation pattern in dB.

## Notes
- The computation time depends on the number of phase shifts and resolution.
- The script automatically generates and executes a Gnuplot script (`plot.gnu`). Ensure Gnuplot is installed for proper visualization.
