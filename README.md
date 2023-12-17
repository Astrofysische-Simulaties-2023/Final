# Contents
### (OUTDATED)

## Getting started

To simply run the C++ code, use NBody_final.cpp. 

## Header Files

Header file |Description 
---|---
[Vec.h](Vec.h) | Contains the 3D vector class.
[Bodyclass.h](Bodyclass.h) | Contains the astronomical body class, with the appropriate methods. We can calculate the acceleration, potential and kinetic energy of a body. In the function calcR, we calculate the closest distance between a specific body and multiple spacecrafts. <br> Adaptive time step: in the function update_dt, we calculate the new time step by using a power law of the normalized acceleration of the bodies. In the function update_dt2, we calculate the new time step by using a max energy error for every simulation time step. We use the forward Euler method for calculating the energy error of a specific time step at a specific time in the simulation. 

[Integrators.h](Integrators.h) | Euler, RK4, Velocity Verlet, Forest Ruth integrators
[ButcherTableau.h](ButcherTableau.h) | Contains the Butcher Tableau class, with different butcher tableaus for a variety of explicit and implicit integrators. Also contains the general functions that transform a butcher tableau into an integrator that updates the position, velocity and acceleration.
[PartionedTableau.h](PartionedTableau.h) | Contains the Partioned Tableau class, with different partioned tableaus for a variety of partioned symplectic integrators. Also contains the general function that transforms a partioned tableau into an integrator that updates the position, velocity and acceleration.
