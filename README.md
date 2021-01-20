## Phormidium pseudopriestleyi MAG repo

This is the directory for the assembling and binning of metagenomic sequencing data of *Phormidium pseudopriestleyi*.

### Directory Structure:
* **Snakefile**: This is a snakeflow of my workflow. It takes a metagenomic sequencing from Lake Fryxell mats and a lab culture of *P. pseudopriestleyi* and assembles them separately, and then as a coassembly. Then anvi'o is used to bin the metagenome assembly.
* **scripts**: These are extra analyses scripts.
* **envs**, **peloton**, and **modules.pel**: These are files necessary for running snakemake on peloton.
