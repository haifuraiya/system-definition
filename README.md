# Phase 4 Space Segment

Top level system definition of the satellite design.

This repository has the following segments:

### Top Level System Modelling

System modelling is performed using the open source [Capella](https://polarsys.org/capella) Model Based Systems Engineering tool to form a dynamic and integrated model of the whole system.

The folder [capella-system-models](capella-system-models) contains the workspace which can be used by the Capella tool.

The outputs of the tool are rendered to navigable HTML, the index of this is found [here](TBD).

### Simulations

System level simulations can be found in the [simulations](simulations) folder.  Although most results are intended to be viewed from the Capella generated HTML.

Current simulations include:
- Solar power generation over specific orbits including radiation effects
- Power budget over orbits including radiation and battery degradation effects
- TID (Total Ionizing Dose) levels expect in SiO 
