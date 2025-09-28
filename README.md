# SPH Fluid + Rigid Body Simulation (C++ / OpenGL)

This repository contains my final project for the IGR202 – Computer Graphics and Virtual Reality class.  
The full project report (PDF) is available in the repo.

The goal was to implement a **2D Smoothed Particle Hydrodynamics (SPH)** fluid coupled with a **rigid body**, rendered in real time using **OpenGL**.  
The project builds on the article *Versatile Rigid–Fluid Coupling for Incompressible SPH* (Akinci et al., 2012) and adapts course material provided by **Pr. Kiwon Um**.

---

## Demo

<img src="final.gif" width="500"/>

---

## Features

- 2D SPH with cubic spline kernel and support radius 2h
- Density and pressure using Tait equation of state
- Simple viscosity term and gravity
- Box wall collisions
- Two-way coupling between fluid and rigid body  
  (forces applied by the fluid on the rigid body are symmetrically applied back to the fluid)
- Rigid body represented by boundary samples
- Uniform grid neighbor search for particles and boundary samples
- Real-time rendering with OpenGL
- Interactive controls:
  - **H**: show help  
  - **P**: play/pause simulation  
  - **G**: toggle grid  
  - **V**: toggle velocity rendering  
  - **S**: save current frame (.tga)  
  - **Q**: quit  

---

## Project Layout

- `main.cpp` — application entry point, render loop, interactive controls
- `SphSolver.hpp` — SPH fluid solver (particles, neighbor grid, density, pressure, viscosity, integration, wall collision, visualization)
- `RigidSolver.hpp` — rigid body solver (state, boundary sampling, force/torque accumulation, two-way coupling)
- `Vector.hpp`, `Vector3.hpp` — small 2D and 3D vector math
- `Matrix3x3.hpp`, `Quaternion.hpp` — minimal rotation math (used for rigid body, coupled back to 2D projection)

---

## Report Highlights

- Implements **two-way fluid–rigid coupling** without particle sticking or penetration issues  
- Handles **boundary contribution** for density and force computations  
- Shows realistic pressure and friction interactions between fluid and solid  
- Performance remains a challenge (neighborhood search is expensive; ~6 min runtime per animation)  
- Possible improvements: neighbor search optimization, merging loops, 3D extension, articulated rigid bodies  

---

## Acknowledgments

- Project supervised by **Pr. Kiwon Um** and **Amal Dev Parakkat**  
- Based on prior course material from Telecom Paris  
- Coupling method inspired by Akinci et al. (2012), *Versatile Rigid–Fluid Coupling for Incompressible SPH*

---
