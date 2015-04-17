**Introduction**

RNA folding pathways are involved in various biological processes, e.g. riboswitches for regulating gene expressions. A challenge to understanding RNA folding reactions is the complex relationships that exist between the structure of the RNA and its folding landscape. The identification of intermediate species that populate conformation changing landscapes and characterization of elements of their structures are the key components to solving the RNA switching problem.

In this study, we implemented the elastic band techniques to model conformational change pathways at the resolution of secondary structure. This method is a chain-of state method in which a series of images is used to describe the path. All the images are connected by spring forces to ensure a similar spacing, in terms of base paring distance, of images along the pathway. The simulation is started from an chain which has to connect reactant and product minimums. A Monte Carlo simulation system is used to sample the pathways in order to identify low free energy change pathways.

**PEB and NEB**

The spring force potential energy between two images E, is determined by the spring constant k and the spring length L, following hook’s law: E = k*(L^2). The length of the spring is defined by base pair difference between two structures of the images. The spring forces potential energy of each image is the sum of the spring forces potential energies of the image and its neighbor images. The total energy of the each image derived from a Monte Carlo step will be evaluated based on the sum of the folding free energy and spring force potential energy. This method is implemented with C++ based on the RNA class in RNAstructure package.
