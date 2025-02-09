# Creator of 2D maps of galaxies from cosmological simulations

A repository for creating 2D spatially-resolved maps of stellar properties of galaxies from cosmological simulations.
Currently, IllustrisTNG and EAGLE cosmological simulations are supported.
These maps can be used for simulation-based inference.

## Background

Cosmological simulations, such as [IllustrisTNG](https://www.tng-project.org/) and [EAGLE](https://eagle.strw.leidenuniv.nl/wordpress/), have transformed our understanding of galaxy formation by producing realistic 3D models of galaxies and large-scale structures. To compare these simulations with observations, we need to create 2D projections of observable properties, mimicking what telescopes capture.

Simulations store properties like stellar mass, luminosity, and metallicity in 3D space. Converting these to 2D maps involves integrating along a chosen line of sight, applying appropriate weighting schemes (e.g., mass or luminosity), and accounting for observational effects such as resolution and orientation. These maps can represent stellar properties, kinematics, or chemical abundances and are key for connecting theory to observation.

This repository provides tools for generating such 2D maps efficiently, enabling direct comparisons with data from galaxy surveys. By bridging simulations and observations, these tools support studies of galaxy structure, evolution, and dynamics.

## Usage

You can install the repo as a dependency in your project:

```bash
pip install -e git+https://github.com/eagel27/cosmological_sims_maps_creator@main#egg=cosmological_sims_maps_creator
```

You can create a provider object in order to run the simulations:

```python
from cosmological_sims_maps_creator.simulations.provider import GalaxySimulationProvider

provider = GalaxySimulationProvider(
        galaxy_type = galaxy_type   # one of "EAGLE", "EAGLE_IllustrisTNG", "IllustrisTNG"
        mode = mode                 # what is mode?
        align = align               # what is align?
        id_start = id_start         # what is id_start?
        id_stop = id_stop           # what is id_stop?
        suffix = suffix             # what is suffix?
)
```

**NOTE:** TODO: Write about assumptions for the data

## Example maps

Maps from IllustrisTNG -- TNG100 (galaxy with id=19 from snapshot=99)
![galaxy_19_99_0](https://github.com/user-attachments/assets/382e954b-9f3f-468d-8bb7-0171046aad8b)

Maps from EAGLE (galaxy with id=19151261 from snapshot=28)
![galaxy_19151261_28_2](https://github.com/user-attachments/assets/6afb574c-fe34-492d-890d-57c5a79bedcf)
