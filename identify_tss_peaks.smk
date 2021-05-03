from snakemake.utils import validate, min_version
import itertools
import pandas as pd

##### set minimum snakemake version #####
min_version("5.20.0")


##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"],sep="\t")

##### target rules #####
def get_input(wildcards):
    if config["analysis"] == "spi":
        robust = "{pdir}/outPooled/tc.spi.merged.ctssMaxCounts{robust}.ctssMaxTpm{tpm}.bed.gz".format(
            pdir=config['directory'],
            robust=config['cutoff']['robust'],
            tpm=config['cutoff']['tpm'])
        permissive = "{pdir}/outPooled/tc.spi.merged.ctssMaxCounts{permissive}.bed.gz".format(
            pdir=config['directory'],
            permissive=config['cutoff']['permissive'])
        return [robust, permissive]
    elif config["analysis"] == "dpi":
        robust = "{pdir}/outPooled/tc.decompose_smoothing_merged.ctssMaxCounts{robust}.ctssMaxTpm{tpm}.bed.gz".format(
            pdir=config['directory'],
            robust=config['cutoff']['robust'],
            tpm=config['cutoff']['tpm'])
        permissive = "{pdir}/outPooled/tc.decompose_smoothing_merged.ctssMaxCounts{permissive}.bed.gz".format(
            pdir=config['directory'],
            permissive=config['cutoff']['permissive'])
        allout = expand("{pdir}/outPooled/tc.long.decompose_smoothing.component{n}_ctss.{strand}.bedGraph.gz",
            pdir=config['directory'],
            strand=['fwd', 'rev'],
            n=range(1, config['split']['n_comp_upper_bound'] + 1))
        allout.append(robust)
        allout.append(permissive)
        return allout
    else:
        raise ValueError("config['analysis'] should be 'spi' or 'dpi'; provided {}".format(config['analysis']))

rule all:
    input:
        get_input
        #"_test/snakemake/outPooled/tc.long/aaaaa"

##### setup report #####
#report: "report/workflow.rst"

##### load rules #####
### prepare bigWig files for individual samples
include: "rules/bigwig.smk"
### pool ctss across samples: total counts, max counts and tpm
include: "rules/pool.smk"
### prepare tag clusters
include: "rules/tagcluster.smk"
### split tag cluster into long and short for decomposition
include: "rules/split_tagclusters.smk"

if config["analysis"] == 'spi':
    ### simple smoothing without decomposition
    include: "rules/smoothing_simple.smk"
    ### thresholding
    include: "rules/thresholding_simple.smk"

if config["analysis"] == 'dpi':
    ### decomposition and smoothing
    include: "rules/smoothing.smk"
    ### thresholding
    include: "rules/thresholding.smk"
