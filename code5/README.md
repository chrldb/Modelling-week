# Antenna Radiation Pattern Simulation

This project simulates the radiation pattern of a phased array antenna by computing the electric field distribution and visualizing the results using Gnuplot.

## Features
- Simulates the electric field propagation for a linear array antenna.
- Computes the received power at different angles.
- Generates radiation pattern plots.
- Uses Gnuplot for visualization.

## Dependencies
To compile and run the program, ensure you have:
- A C++ compiler (GCC or Clang recommended)
- Gnuplot installed for plotting

## Compilation & Execution

### Compile the program:
```sh
 g++ -o antenna_simulation antenna_simulation.cpp -std=c++11 -O2
```

### Run the simulation:
```sh
 ./antenna_simulation
```

## Output
- The program generates a data file `data_cartesian.dat` containing computed values.
- A Gnuplot script `cartesian_plot.gnu` is created.
- The final visualization is stored in `cartesian_plot.png`.

## How It Works
### 1. Antenna Model
- Defines a linear array antenna with a configurable number of emitters.
- Uses a phase shift `delta_phi` to steer the main beam.

### 2. Computation
- Calculates the electric field at each observation point.
- Computes the received power at different angles.
- Normalizes the power and converts it to dB scale.

### 3. Visualization
- Outputs the computed data to a file.
- Uses Gnuplot to generate the radiation pattern graph.

## Example Configuration
Default parameters in the code:
- **Operating frequency:** 1.2 GHz
- **Number of emitters:** 40
- **Delta phase shift:** π/4
- **Scan radius:** 1000m

## Notes
- The script automatically opens the generated image (`cartesian_plot.png`).
- Ensure Gnuplot is correctly installed and accessible via command line.
