# Creator of 2D maps of galaxies from cosmological simulations 
A repository for creating 2D spatially-resolved maps of stellar properties of galaxies from cosmological simulations. 
Currently, IllustrisTNG and EAGLE cosmological simulations are supported. 
These maps can be used for simulation-based inference.

## Background
Cosmological simulations, such as [IllustrisTNG](https://www.tng-project.org/) and [EAGLE](https://eagle.strw.leidenuniv.nl/wordpress/), have transformed our understanding of galaxy formation by producing realistic 3D models of galaxies and large-scale structures. To compare these simulations with observations, we need to create 2D projections of observable properties, mimicking what telescopes capture.

Simulations store properties like stellar mass, luminosity, and metallicity in 3D space. Converting these to 2D maps involves integrating along a chosen line of sight, applying appropriate weighting schemes (e.g., mass or luminosity), and accounting for observational effects such as resolution and orientation. These maps can represent stellar properties, kinematics, or chemical abundances and are key for connecting theory to observation.

This repository provides tools for generating such 2D maps efficiently, enabling direct comparisons with data from galaxy surveys. By bridging simulations and observations, these tools support studies of galaxy structure, evolution, and dynamics.

## Usage





