# Chaotic pendulum Simulation (Finite Difference Method)
This project, carried out as a group assignment for the **LPHYS1303** course, implements a numerical simulation of the **double pendulum**, a nonlinear dynamic system known for its **chaotic behaviour**. The equations of motion are integrated using the **fourth-order Runge-Kutta method (RK4)**.
## Description
The double pendulum consists of two masses suspended from each other. Its equations of motion are obtained using **Lagrangian formalism** and exhibit **non-linear coupling**, making the system very sensitive to initial conditions.

The objective is to:
- Numerically solve these equations with **Runge-Kutta 4**
- To qualitatively analyse different trajectories
- To study **sensitivity to initial conditions** and **chaos** through the Lyapunov exponent
## Use
Run the code to launch the simulation: [Pendule_Double.py](Pendule_Double.py)
## Results
- **Simulation of trajectories:** analysis of the movement of masses for different initial conditions.
- **Evidence of chaos:** visualisation of the evolution of angles and velocities.
- **Quantitative analysis:** study of the **Lyapunov** exponent to characterise the divergence of trajectories.
