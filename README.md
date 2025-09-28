# SPH Fluid + Rigid Body Simulation (C++ / OpenGL)

This repository contains my final project for the IGR202-Computer Graphics and VR class.

I implemented a small real time solver to create a 2D Smoothed Particle Hydrodynamics (SPH) fluid coupled with a rigid body, rendered with OpenGL. 

Note : Parts of the code were written by Pr. Um Kiwon which I adapted for this project.


##Demo
<img src="final.gif" width="500"/>


## Features

- 2D SPH with cubic spline kernel and support radius 2h
- Density and pressure using Tait equation of state, simple viscosity term
- Gravity and box wall collisions
- Rigid body boundary represented by sample points, two way coupling with the fluid
- Uniform grid neighbor search for particles and boundary samples
- Immediate mode style OpenGL rendering for simplicity
- Cross platform build with CMake

## Project layout
- `SphSolver.hpp` - 2D SPH solver: particles, neighbor grid, density, pressure, viscosity, integration, wall collision, colors and velocity lines for visualization.
- `RigidSolver.hpp` - Rigid body state and boundary sampling. Builds world space positions for boundary points and the map cell -> indices for neighbor queries. Accumulates forces and torques from the fluid.
- `Vector.hpp`, `Vector3.hpp` - small 2D and 3D vector math.
- `Matrix3x3.hpp`, `Quaternion.hpp` - minimal 3D rotation math used by the rigid body. Final coupling to SPH uses 2D projections.
- `main.cpp` - window and render loop. 

