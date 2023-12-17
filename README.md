# Contents


To run the C++ code, use NBody_final.cpp. <br>
To use adaptive timestep, uncomment line 128 or 129 (for the first, or the second method).

## Header Files

Header file |Description 
---|---
[Vec.h](Vec.h) | 3D vector class.
[Bodyclass.h](Bodyclass.h) | Contains the astronomical body class, with the appropriate methods (acceleration, energy). <br> Adaptive time step: two methods are implemented. The first using a power of the acceleration. The second trying to minimize the energy difference, using a Forward Euler method.
[Integrators.h](Integrators.h) | Euler, RK4, Velocity Verlet, Forest Ruth integrators
[ButcherTableau.h](ButcherTableau.h) | Contains the Butcher Tableau class, with different butcher tableaus for several explicit and implicit integrators. Also contains the general functions that transform a butcher tableau into an integrator. Implicit RK3 and RK5 are implemented in the final program.
[Choose_integrator.h](Choose_integrator.h) | Methods to choose the integrator at user input.

