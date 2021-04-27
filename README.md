Decomposition-based peak identification (DPI)
=============================================

This repository is a re-implementation of DPI, originally written by
(hkawaji)[https://github.com/hkawaji/dpi1].

DPI identifies a set of reference peaks in in genome-wide TSS (transcription
starting site) profiles obtained from diverse range of biological states. In
particular, this is developed for the FANTOM5 project (Nature 507, 462-467, 2014
and its related papers) that produced CAGE profiles in more than one thousand of
biological states in human and mouse, and the analysis results have been
available at http://fantom.gsc.riken.jp/5/ as well as the article above.

The original implementation uses a combination of bash scripts, RAKE and
and integrated compute cluster commands (via qsub). Here, we use the
snakemake workflow manager with integrated conda environments to implement DPI.
It can be run on a local computer or making use of a high performance compute
cluster.

Requirements
------------

  - bash (https://www.gnu.org/software/bash/)
  - ruby (https://www.ruby-lang.org)
  - R (http://cran.r-project.org/)
  - R/fastICA package (http://cran.r-project.org/web/packages/fastICA/)
  - R/tidyr package (https://cran.r-project.org/web/packages/tidyr/)
  - command line tools to operate bigWig files in http://hgdownload.cse.ucsc.edu/admin/
  - bedTools (https://code.google.com/p/bedtools/)
  - Unix/Linux with Grid Engine (developed with UGE)

Installation
------------

    % git clone https://github.com/hkawaji/dpi1.git

How to run
-----------

    % cd ${INSTALLED_DIR}
    % ./identify_tss_peaks.sh [options]

Follow the usage described in the message.




Author
------
Hannah V Meyer

Reference
---------
* A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. Nature 507, 462-467, 2014


