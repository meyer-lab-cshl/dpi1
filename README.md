Decomposition-based peak identification (DPI)
=============================================

This repository is a re-implementation of DPI, originally written by
[hkawaji](https://github.com/hkawaji/dpi1).

DPI identifies a set of reference peaks in genome-wide TSS (transcription
starting site) profiles. In particular, this was developed for the FANTOM5
project (Nature 507, 462-467, 2014 and its related papers) that produced CAGE
profiles in more than one thousand of biological states in human and mouse, and
the analysis results have been available at
http://fantom.gsc.riken.jp/5/ as well as the article above.

The original implementation uses a combination of bash scripts, RAKE and
and integrated compute cluster commands (via qsub). Here, we use the
snakemake workflow manager with integrated conda environments to implement DPI.
It can be run on a local computer or making use of a high performance compute
cluster.

Requirements
------------
Users need to have access to a unix system and should install snakemake.
Please follow these [detailed install
instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

This workflow uses conda as a package manager. For each rule where external
software is required, an conda environment is provided which will be
automatically installed when running the workflow for the first time.

In addition there are some tools not accessible by conda. Please install
the following prior to running the workflow and make sure they are located
in your PATH
  - [bedTools](https://code.google.com/p/bedtools/): mergeBed, intersectBed

Installation
------------

    git clone https://github.com/meyer-lab-cshl/dpi1

How to run
-----------
1. Move to the DPI directory:
  
  ```cd ${INSTALLED_DIR}```

2. Adapt the config/config.yaml file with your system specific and analysis
parameters. In particular, these parameters have to be set:

* Type of analysis: spi/dpi (simple peak calling without ICA-based decomposition
or peak calling with ICA-based decomposition as in orignal paper)
    * analysis: dpi

* Files and directories: use full file paths, do NOT use tilde expansion
    
    * directory: path/to/data/directory
    * genome: path/to/organism_genome_file_for_bedtools (containing chromosome
    * name and length, for example see _test/mouse.mm10.genome)
    * samples:
        - sample_1
        - sample_2

(This assumes sample files having the suffix .ctss.bed.gz, eg
sample_1.ctss.bed.gz)

The remaining parameters are the default parameters from the original
implementation.

3. Activate your snakemake environment
4. Run (local machine, using 1 core):

    ```snakemake -s identify_tss_peaks.smk --cores 1 --use-conda```


Author
------
Hannah V Meyer

Reference
---------
* A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. Nature 507, 462-467, 2014


