# Interacting Gas Simulation

This repository contains Python code for simulating an interacting gas system using **Hamiltonian dynamics** and a **semi-implicit (symplectic) Euler integration scheme**. The project is intended as a numerical and educational tool for studying many-particle dynamics, energy-conserving time integration, and emergent behavior in classical systems.


---

## Overview

The simulation models a gas composed of interacting particles evolving according to a Hamiltonian formulation. The equations of motion are derived from a specified Hamiltonian (typically consisting of kinetic and potential energy terms), and the system is evolved forward in time using a semi-implicit Euler method, which offers improved stability and qualitative energy behavior compared to the explicit Euler method.

Key features include:
- Classical Hamiltonian dynamics for many-body systems
- Pairwise interaction potentials between particles
- Semi-implicit (symplectic) Euler time integration
- Configurable number of particles, dimension, time step, and interaction parameters
- Basic diagnostics such as energy evolution and particle trajectories

---

## Mathematical Model

### Hamiltonian Dynamics

The system is governed by a Hamiltonian of the form

\[ H(q, p) = T(p) + V(q) \]

where:
- \( q \) denotes particle positions,
- \( p \) denotes particle momenta,
- \( T(p) \) is the kinetic energy,
- \( V(q) \) is the potential energy.

The potential energy is
The equations of motion are given by Hamilton’s equations:

\[
\dot{q} = \frac{\partial H}{\partial p}, \qquad
\dot{p} = -\frac{\partial H}{\partial q}.
\]

---

### Semi-Implicit Euler Method

Time integration is performed using the semi-implicit (symplectic) Euler scheme:

\[
\begin{aligned}
 p_{n+1} &= p_n - \Delta t \, \nabla_q V(q_n), \\
 q_{n+1} &= q_n + \Delta t \, \nabla_p T(p_{n+1}).
\end{aligned}
\]

This method is symplectic, making it well-suited for Hamiltonian systems where long-term qualitative behavior (such as approximate energy conservation) is important.

---

## Repository Structure

```
.
├── src/            # Core simulation code
├── examples/       # Example scripts and configurations
├── utils/          # Helper functions (initial conditions, plotting, etc.)
├── README.md       # Project documentation
└── requirements.txt
```

*(Directory names may vary depending on the version of the project.)*

---

## Requirements

- Python 3.8 or newer
- NumPy
- Matplotlib (for visualization)

Install dependencies with:

```bash
pip install -r requirements.txt
```

---

## Usage

A typical workflow consists of:
1. Defining the number of particles and interaction parameters
2. Initializing positions and momenta
3. Running the time integration loop
4. Analyzing or visualizing the results

Example:
The notebook: testing\_notebook.ipynb contains an example how to use the code
---

## Output and Visualization

The simulation can produce:
- Particle trajectories
- Phase-space plots
- Energy vs. time plots (to assess numerical stability)

Visualization scripts are provided in the `examples/` or `utils/` directories.

---

## Limitations and Notes

- The semi-implicit Euler method is conditionally stable and still requires sufficiently small time steps.
- The simulation is classical and does not include quantum effects.
- Performance may degrade for very large particle counts without further optimization (e.g., neighbor lists or parallelization).

---

## Intended Audience

This repository is suitable for:
- Students learning computational physics or molecular dynamics
- Researchers prototyping Hamiltonian systems
- Anyone interested in symplectic integrators and many-body dynamics

---

## License

Specify the license under which this project is distributed (e.g., MIT, BSD, GPL). If no license is provided, all rights are reserved by default.

---

