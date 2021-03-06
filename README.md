# mccloud
![Clumpy](/clumpyPhase.PNG)
Reverse Mote Carlo approach to modeling ray propagation through a clumpy molecular cloud.

## Introduction
Molecular clouds are interstellar clouds of size and density capable of forming molecules. Observations show the molecules are not evenly distributed throughout the cloud’s structure. The high density "clumps" are the beginning of star formation. The chemical reaction rates driving the stellar development are dependent on the ambient radiation field. Therefore, an understanding of radiation within molecular clouds is important to understanding stellar evolution. This project will investigate the radiation field by first modeling a structured molecular cloud and then preforming a reverse Monte-Carlo procedure to quantify the propagation of radiation from changing views and density contrasts.

## Implementation
The following code was written for a class project during the fall 2017 semester at UNC-Chapel Hill (PHYS 358: Physical Modeling). The implementation is based on the paper *Dust heating in turbulent molecular cloud models* (Bethell et al. 2004).

For a description of the model and discussion of the results, please see the [final report](https://github.com/Ryan-A-Armstrong/mccloud/blob/master/FINALMCpaper.pdf) and [class presentation](https://github.com/Ryan-A-Armstrong/mccloud/blob/master/MCCloudPresentation.pdf).

## Running the code
Instructions for running the code are provided in comments at the beginning of the two .py files.

mccloud_transfer.py creates the model and runs the simulation. 

analysis.py allows for visualizaiton of the results of the simulation. 
