# Antenna Simulation Project

## Overview
This repository contains multiple C++ programs related to antenna modeling and radar beamforming. The project aims to simulate different aspects of phased-array antennas, including beamforming, power calculations, phase shifting, and advanced windowing techniques.

## Structure of the Repository
The repository includes the following main programs:

### 1. **Gain Calculation Based on Theta**
   - Computes the gain of the antenna as a function of the angle \( \theta \) for a given radius.
   - Helps analyze how the gain changes with the observation angle.
   - Outputs data that can be plotted using Gnuplot.

### 2. **Far-Field Power Calculation**
   - Uses far-field approximation to calculate the received power at a given distance.
   - Implements the classic radar equation to analyze power loss over long distances.
   - Helps understand the impact of distance on power reception.

### 3. **Phase Shift Modeling for Scanning a Given Range**
   - Determines the number of phase shifts required to scan a given angular range (e.g., \(-\pi/4\) to \(+\pi/4\)).
   - Ensures proper beam coverage with overlapping lobes at \(-3\) dB.
   - Simulates scanning time based on round-trip travel time of radar signals.

### 4. **Antenna Weighting with Blackman, Chebyshev, and Taylor Windows**
   - Implements different windowing techniques to shape the antenna radiation pattern.
   - Reduces sidelobes using:
     - **Blackman** window
     - **Chebyshev** window (with a given ripple factor)
     - **Taylor** window (optimized for side lobe level control)
   - Useful for improving target detection by suppressing unwanted sidelobes.

## Usage
Each program is compiled and executed separately. Below is an example compilation command for a typical program:

```bash
 g++ -o antenna_simulation antenna_simulation.cpp -O2 -std=c++11 -lm
 ./antenna_simulation <parameters>
```

### Example Run Commands
#### **Gain Calculation**
```bash
./antenna_gain 40 1.2e9
```
#### **Far-Field Power Calculation**
```bash
./far_field_power 40 1.2e9 0.1 0.5
```
#### **Phase Shift Simulation**
```bash
./phase_shift_simulation 40 1.2e9 -pi/4 pi/4
```
#### **Windowed Beamforming**
```bash
./beamforming_window 40 1.2e9 chebyshev -30.0
```

## Visualization
Each program generates output files and a .png picture using **Gnuplot**. 
