#!/usr/bin/env snakemake
# File: association_analysis.smk
# Author: Zhenhua Zhang
# Created: Jul 11, 2023
# Updated: Jul 27, 2023

"""
This is a snakemake pipeline for the elite controller rare-variant association analysis.

1. Required tools
  - regenie
  - plink2.0

2. Primer files
# 2.1 Prepare variant set of gene
python3 scripts/py3/prepare_variant_sets.py \
  -p 4 \
  -b outputs/associations/setbased/whole_genome_genes.bed \
  -o outputs/associations/setbased/whole_genome_genes \
  outputs/variants/gsa/annotation/per_chr/2000HIV.all.ann.vcf.gz
"""

import csv
import pathlib

import yaml
from benedict import benedict
# from gwastk import PlinkCmd, RegenieCmd


# Configuration using YAML. Please pass the custom configuration file by `--config myconfig=PATH/TO/CONFIG/FILE`
class AsPath:
  def __init__(self, *values):
    self.values = values
  def parse(self):
    return pathlib.Path("/".join([str(x) for x in self.values])).expanduser()


def as_path(loader, node):
  """A custom YAML constructor for Path"""
  if isinstance(node, yaml.ScalarNode):
    values = loader.construct_scalar(node)
    return AsPath(values).parse()
  elif isinstance(node, yaml.SequenceNode):
    values = loader.construct_sequence(node)
    return AsPath(*values).parse()
  raise NotImplementedError("Unsupported node type: {}".format(type(node)))


# Obtain categorical covariates
def get_cat_var_list(wildcards, output):
  return ",".join(bconfig.get("cat_covars/default", []) + bconfig.get(f"cat_covars/{wildcards.per_partition}", []))


yaml.add_constructor("!AsPath", as_path)

bconfig = None
with open(config["myconfigfile"], "r") as f:
  config = yaml.load(f, Loader=yaml.FullLoader)
  bconfig = benedict(config, keypath_separator = "/")


#
## Resources
#
# Project directory
project_dir = bconfig.get("project_dir")
project_token = bconfig.get("project_token")

# Singularity containers
singularity_image = bconfig.get("singularity_image")

# Genotype, phenotype
genotype_dir = bconfig.get("genotype_dir")
phenotype_dir = bconfig.get("phenotype_dir")

# Covariables. Covariables, interaction covaraibles
covariate_file = bconfig.get("covariate_file")

# Interaction with E/G, i.e., GxE or GxG
interaction_with = bconfig.get("interaction_with")
if interaction_with is None:
  regenie_int_opt = ""
  plink_split_cat_pheno = ""
  plink_int_opt = ""
else:
  regenie_int_opt = "--interaction " + interaction_with
  plink_int_opt = "interaction"# --parameters " + bconfig.get("params/plink/parameters")
  plink_cat_var_list = bconfig.get("cat_covars/default", []) + bconfig.get("cat_covars/plink", [])
  if plink_cat_var_list:
    plink_split_cat_pheno = "--split-cat-pheno " + " ".join(plink_cat_var_list)
  else:
    plink_split_cat_pheno = ""


# Using variants as covariables. Only works for regenie.
condition_variant = bconfig.get("condition_variant")
if condition_variant is not None:
  pass # TODO

# Regression model, logistic (GWAS) or linear (QTL)
regression_model = bconfig.get("model")

# Groups of individuals. all, discovery, and validation
partition_list = bconfig.get("partition")

# Which phenotype to work on.
phenotype_list = bconfig.get("phenotype")

# Global configures
plink_vif = bconfig.get("params/plink/vif")
plink_maf = bconfig.get("params/default/maf")
plink_mac = bconfig.get("params/default/mac")
plink_out_suffix = "logistic.hybrid" if regression_model == "logistic" else "linear"
plink_firth_fb = "firth-fallback" if regression_model == "logistic" else ""

regenie_aaf_bins = ",".join([str(x) for x in bconfig.get("params/regenie/aaf_bins")])
regenie_mhc_bed_file = bconfig.get("params/regenie/mhc_regions") # To be excluded, major histocompatibility complex region
regenie_nlc_bed_file = bconfig.get("params/regenie/nlc_regions") # To be included, NON low complexity regions
regenie_max_cat_levels = bconfig.get("params/regenie/max_cat_levels")
regenie_min_mac = bconfig.get("params/default/mac")
regenie_min_maf = bconfig.get("params/default/maf")
regenie_model_flag = "--bt --af-cc" if regression_model == "logistic" else "--qt"


# Outputs
final_output = []
output_dir = bconfig.get("output_dir")

plink_gwas_out_file = expand(output_dir/"gwas/summary_statistic/{per_pheno}/{per_partition}"/("gwas_sumstat.{per_pheno}.glm." + plink_out_suffix + ".gz"), per_pheno=phenotype_list, per_partition=partition_list)
plink_xqtl_out_file = expand(output_dir/"xqtl/summary_statistic/{per_pheno}/{per_partition}"/("gwas_sumstat.{per_pheno}.glm." + plink_out_suffix + ".gz"), per_pheno=phenotype_list, per_partition=partition_list)
plinK_setbased_out_file = expand(output_dir/"setbased/summary_statistic/{per_pheno}/{per_partition}"/("setbased_sumstat.{per_pheno}.glm." + plink_out_suffix + ".gz"), per_pheno=phenotype_list, per_partition=partition_list)
regenie_gwas_out_file = expand(output_dir/"gwas/summary_statistic/{per_pheno}/{per_partition}/gwas_sumstat_{per_pheno}.regenie.gz", per_pheno=phenotype_list, per_partition=partition_list)
regenie_setbased_out_file = expand(output_dir/"setbased/summary_statistic/{per_pheno}/{per_partition}/setbased_sumstat_{per_pheno}.regenie.gz", per_pheno=phenotype_list, per_partition=partition_list)
regenie_gwas_simple_out_file = expand(output_dir/"gwas/summary_statistic/{per_pheno}/{per_partition}/gwas_sumstat.simple_{per_pheno}.regenie.gz", per_pheno=phenotype_list, per_partition=partition_list)
qtltools_xqtl_out_file = expand(output_dir/"xqtl/summary_statistic/{per_pheno}/{per_partition}/xqtl_sumstat.{per_pheno}.qtltools.gz", per_pheno=phenotype_list, per_partition=partition_list)

run_mode_plink, run_mode_regenie, run_mode_qtltools = bconfig.get("run_mode/plink", []), bconfig.get("run_mode/regenie", []), bconfig.get("run_mode/qtltools", [])
if run_mode_plink is not None and "gwas" in run_mode_plink: final_output += plink_gwas_out_file
if run_mode_plink is not None and "setbased" in run_mode_plink: final_output += plink_setbased_out_file
if run_mode_plink is not None and "xqtl" in run_mode_plink: final_output += plink_xqtl_out_file
if run_mode_regenie is not None and "gwas_simple" in run_mode_regenie: final_output += regenie_gwas_simple_out_file
if run_mode_regenie is not None and "gwas" in run_mode_regenie: final_output += regenie_gwas_out_file
if run_mode_regenie is not None and "setbased" in run_mode_regenie: final_output += regenie_setbased_out_file
if run_mode_qtltools is not None and "xqtl" in run_mode_qtltools: final_output += qtltools_xqtl_out_file


# Wildcard constraints
wildcard_constraints:
  per_pheno = "\w+"


#
## Rules
#
rule all:
  input: final_output


rule s01_perform_gwas_plink:
  input:
    pheno_file = phenotype_dir/(project_token + ".{per_pheno}.txt"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
    covar_file = covariate_file,
    bed_file = genotype_dir/"plink"/(project_token + ".all.bed"),
    bim_file = genotype_dir/"plink"/(project_token + ".all.bim"),
    fam_file = genotype_dir/"plink"/(project_token + ".all.fam"),
  output:
    output_dir/"gwas/summary_statistic/{per_pheno}/{per_partition}/gwas_sumstat.{per_pheno}.glm.{plink_out_suffix}.gz",
  params:
    in_prefix = lambda wildcards, input: pathlib.Path(input.bed_file).with_suffix(""),
    out_prefix = lambda wildcards, output: output[0].replace(".glm", "").replace(f".{plink_out_suffix}.gz", ""),
    out_cols = "chrom,pos,ref,alt,provref,omitted,a1count,totallele,a1countcc,totallelecc,gcountcc,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p,err"
  resources:
    cpus_per_task = 8, mem = "8G", time = "24:00:00"
  shell:
    r"""
    set +eo

    mkdir -p $(dirname {params.out_prefix})
    apptainer exec -B /vol {singularity_image} plink2.0 \
      --glm hide-covar allow-no-covars cols={params.out_cols} {plink_firth_fb} {plink_int_opt} \
      --1 --no-sex --covar-variance-standardize --write-covar --write-samples \
      --vif {plink_vif} \
      --maf {plink_maf} \
      --mac {plink_mac} \
      --bfile {params.in_prefix} \
      --keep {input.kept_donor_file} \
      --pheno {input.pheno_file} \
      --pheno-name {wildcards.per_pheno} \
      --covar {input.covar_file} \
      --threads {resources.cpus_per_task} \
      --out {params.out_prefix}

    gzip -c {params.out_prefix}.{wildcards.per_pheno}.glm.{plink_out_suffix} >| {output}
    rm -f {params.out_prefix}.{wildcards.per_pheno}.glm.{plink_out_suffix}
    """


rule s01_perform_xqtl_plink:
  input:
    pheno_file = phenotype_dir/(project_token + ".{per_pheno}.txt"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
    covar_file = covariate_file,
    bed_file = genotype_dir/"plink"/(project_token + ".all.bed"),
    bim_file = genotype_dir/"plink"/(project_token + ".all.bim"),
    fam_file = genotype_dir/"plink"/(project_token + ".all.fam"),
  output:
    output_dir/"xqtl/summary_statistic/{per_pheno}/{per_partition}/gwas_sumstat.{per_pheno}.glm.{plink_out_suffix}.gz",
  params:
    in_prefix = lambda wildcards, input: pathlib.Path(input.bed_file).with_suffix(""),
    out_prefix = lambda wildcards, output: output[0].replace(".glm", "").replace(f".{plink_out_suffix}.gz", ""),
    out_cols = "chrom,pos,ref,alt,provref,omitted,a1count,totallele,a1countcc,totallelecc,gcountcc,a1freq,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p,err"
  resources:
    cpus_per_task = 8, mem = "8G", time = "24:00:00"
  shell:
    r"""
    set +eo

    mkdir -p $(dirname {params.out_prefix})
    apptainer exec -B /vol {singularity_image} plink2.0 \
      --glm hide-covar allow-no-covars cols={params.out_cols} {plink_firth_fb} {plink_int_opt} \
      --no-sex --covar-variance-standardize --write-covar --write-samples \
      --vif {plink_vif} \
      --maf {plink_maf} \
      --mac {plink_mac} \
      --bfile {params.in_prefix} \
      --keep {input.kept_donor_file} \
      --pheno {input.pheno_file} \
      --pheno-name {wildcards.per_pheno} \
      --covar {input.covar_file} \
      --threads {resources.cpus_per_task} \
      --out {params.out_prefix}

    gzip -c {params.out_prefix}.{wildcards.per_pheno}.glm.{plink_out_suffix} >| {output}
    rm -f {params.out_prefix}.{wildcards.per_pheno}.glm.{plink_out_suffix}
    """


rule s01_perform_gwas_regenie_simple:
  input:
    pheno_file = phenotype_dir/(project_token + ".{per_pheno}.txt"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
    covar_file = covariate_file,
    bgen_file = genotype_dir/"bgen"/(project_token + ".all.bgen"),
    bgi_file = genotype_dir/"bgen"/(project_token + ".all.bgen.bgi"),
    sample_file = genotype_dir/"bgen"/(project_token + ".all.sample"),
  output:
    output_dir/"gwas/summary_statistic/{per_pheno}/{per_partition}/gwas_sumstat.simple_{per_pheno}.regenie.gz",
  params:
    b_size = 200,
    out_prefix = lambda wildcards, output: output[0].replace(f"_{wildcards.per_pheno}.regenie.gz", ""),
    cat_covar_list = get_cat_var_list
  resources:
    cpus_per_task = 8, mem = "32G", time = "8:00:00"
  shell:
    r"""
    set +eo

    mkdir -p $(dirname {params.out_prefix})
    apptainer exec -B /vol {singularity_image} regenie \
      --step 2 \
      --ref-first \
      --keep {input.kept_donor_file} \
      --bgen {input.bgen_file} \
      --bgi {input.bgi_file} \
      --sample {input.sample_file} \
      --phenoFile {input.pheno_file} \
      --phenoColList {wildcards.per_pheno} \
      --covarFile {input.covar_file} \
      --catCovarList {params.cat_covar_list} \
      --maxCatLevels {regenie_max_cat_levels} \
      --threads {resources.cpus_per_task} \
      --minMAC {regenie_min_mac} \
      --bsize {params.b_size} \
      {regenie_int_opt} \
      {regenie_model_flag} --gz --af-cc --write-samples --ignore-pred \
      --out {params.out_prefix}
    """


rule s01_select_variants: # Select genetic variants for whole genome regression model
  input:
    bed_file = genotype_dir/"plink"/(project_token + ".all.bed"),
    bim_file = genotype_dir/"plink"/(project_token + ".all.bim"),
    fam_file = genotype_dir/"plink"/(project_token + ".all.fam"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
  output:
    incl_marker_file = output_dir/"model/{per_pheno}/{per_partition}/step1_variants.prune.in",
    excl_marker_file = output_dir/"model/{per_pheno}/{per_partition}/step1_variants.prune.out",
  params:
    in_prefix = lambda wildcards, input: pathlib.Path(input.bed_file).with_suffix(""),
    out_prefix = lambda wildcards, output: output.incl_marker_file.replace(".prune.in", ""),
  resources:
    cpus_per_task = 4, mem = "8G", time = "1:59:00"
  shell:
    r"""
    set +eo

    apptainer exec -B /vol {singularity_image} plink1.9 \
      --bfile {params.in_prefix} \
      --maf 0.1 \
      --hwe 1e-4 \
      --geno 0.1 \
      --indep-pairwise 1000 100 0.8 \
      --keep {input.kept_donor_file} \
      --exclude range {regenie_mhc_bed_file} \
      --extract range {regenie_nlc_bed_file} \
      --threads {resources.cpus_per_task} \
      --out {params.out_prefix}
    """


rule s02_fit_wgr_model: # Step 1. The whole genome regression model is fit to the traits
  input:
    pheno_file = phenotype_dir/(project_token + ".{per_pheno}.txt"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
    covar_file = covariate_file,
    bed_file = genotype_dir/"plink"/(project_token + ".all.bed"),
    bim_file = genotype_dir/"plink"/(project_token + ".all.bim"),
    fam_file = genotype_dir/"plink"/(project_token + ".all.fam"),
    marker_file = output_dir/"model/{per_pheno}/{per_partition}/step1_variants.prune.in",
  output:
    output_dir/"model/{per_pheno}/{per_partition}/step1_model.{per_pheno}_pred.list",
  params:
    b_size = 1000,
    in_prefix = lambda wildcards, input: pathlib.Path(input.bed_file).with_suffix(""),
    out_prefix = lambda wildcards, output: output[0].replace("_pred.list", ""),
    cat_covar_list = get_cat_var_list,
    model_flag = regenie_model_flag if "--qt" in regenie_model_flag else regenie_model_flag + " --spa"
  resources:
    cpus_per_task = 8, mem = "32G", time = "2:59:00"
  shell:
    r"""
    set +eo

    apptainer exec -B /vol {singularity_image} regenie \
      --step 1 \
      --keep {input.kept_donor_file} \
      --bed {params.in_prefix} \
      --phenoFile {input.pheno_file} \
      --phenoColList {wildcards.per_pheno} \
      --covarFile {input.covar_file} \
      --catCovarList {params.cat_covar_list} \
      --maxCatLevels {regenie_max_cat_levels} \
      --extract {input.marker_file} \
      --threads {resources.cpus_per_task} \
      --minMAC {regenie_min_mac} \
      --bsize {params.b_size} \
      {params.model_flag} --loocv \
      --out {params.out_prefix}
    """


rule s03_perform_gwas_regenie: # Step 2. Perform Firth association using all SNPs
  input:
    model_list = output_dir/"model/{per_pheno}/{per_partition}/step1_model.{per_pheno}_pred.list",
    pheno_file = phenotype_dir/(project_token + ".{per_pheno}.txt"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
    covar_file = covariate_file,
    bgen_file = genotype_dir/"bgen"/(project_token + ".all.bgen"),
    bgi_file = genotype_dir/"bgen"/(project_token + ".all.bgen.bgi"),
    sample_file = genotype_dir/"bgen"/(project_token + ".all.sample"),
  output:
    output_dir/"gwas/summary_statistic/{per_pheno}/{per_partition}/gwas_sumstat_{per_pheno}.regenie.gz",
  params:
    b_size = 200,
    out_prefix = lambda wildcards, output: output[0].replace(f"_{wildcards.per_pheno}.regenie.gz", ""),
    cat_covar_list = get_cat_var_list,
    model_flag = regenie_model_flag if "--qt" in regenie_model_flag else regenie_model_flag + " --spa"
  resources:
    cpus_per_task = 8, mem = "32G", time = "8:00:00"
  shell:
    r"""
    set +eo

    mkdir -p $(dirname {params.out_prefix})
    apptainer exec -B /vol {singularity_image} regenie \
      --step 2 \
      --ref-first \
      --keep {input.kept_donor_file} \
      --bgen {input.bgen_file} \
      --bgi {input.bgi_file} \
      --sample {input.sample_file} \
      --phenoFile {input.pheno_file} \
      --phenoColList {wildcards.per_pheno} \
      --covarFile {input.covar_file} \
      --catCovarList {params.cat_covar_list} \
      --maxCatLevels {regenie_max_cat_levels} \
      --pred {input.model_list} \
      --threads {resources.cpus_per_task} \
      --minMAC {regenie_min_mac} \
      --bsize {params.b_size} \
      {regenie_int_opt} \
      {params.model_flag} --gz --af-cc --write-samples \
      --out {params.out_prefix}
    """


rule s03_perform_setbased_regenie: # Step 2. Perform gene (or any unit) based test using genomic units compiled with rare-variants
  input:
    model_list = output_dir/"model/{per_pheno}/{per_partition}/step1_model.{per_pheno}_pred.list",
    pheno_file = phenotype_dir/(project_token + ".{per_pheno}.txt"),
    kept_donor_file = phenotype_dir/(project_token + ".{per_pheno}.{per_partition}.txt"),
    covar_file = covariate_file,
    bgen_file = genotype_dir/"bgen"/(project_token + ".all.bgen"),
    bgi_file = genotype_dir/"bgen"/(project_token + ".all.bgen.bgi"),
    sample_file = genotype_dir/"bgen"/(project_token + ".all.sample"),
    var_annotation_file = output_dir/"setbased/whole_genome_genes.var_annotation.txt",
    var_setlist_file = output_dir/"setbased/whole_genome_genes.var_set_list.txt",
    var_maks_file = output_dir/"setbased/whole_genome_genes.mask.txt",
  output:
    output_dir/"setbased/summary_statistic/{per_pheno}/{per_partition}/setbased_sumstat_{per_pheno}.regenie.gz"
  params:
    b_size = 200,
    out_prefix = lambda wildcards, output: output[0].replace(f"_{wildcards.per_pheno}.regenie.gz", ""),
    cat_covar_list = get_cat_var_list,
    model_flag = regenie_model_flag
  resources:
    cpus_per_task = 8, mem = "32G", time = "8:00:00"
  shell:
    r"""
    set +eo

    mkdir -p $(dirname {params.out_prefix})
    apptainer exec -B /vol {singularity_image} regenie \
      --step 2 \
      --ref-first \
      --keep {input.kept_donor_file} \
      --bgen {input.bgen_file} \
      --bgi {input.bgi_file} \
      --sample {input.sample_file} \
      --phenoFile {input.pheno_file} \
      --phenoColList {wildcards.per_pheno} \
      --covarFile {input.covar_file} \
      --catCovarList {params.cat_covar_list} \
      --maxCatLevels {regenie_max_cat_levels} \
      --pred {input.model_list} \
      --anno-file {input.var_annotation_file} \
      --set-list {input.var_setlist_file} \
      --mask-def {input.var_maks_file} \
      --aaf-bins {regenie_aaf_bins} \
      --threads {resources.cpus_per_task} \
      --bsize {params.b_size} \
      --minMAC {regenie_min_mac} \
      --minINFO 0.9 \
      {regenie_int_opt} {params.model_flag} --gz --write-samples --write-mask --write-mask-snplist \
      --out {params.out_prefix}
    """


rule s03_perform_xqtl_qtltools:
  input:
    vcf_file = output_dir/"",
    vcf_tbi_file = output_dir/"",
    bed_file = output_dir/"",
    covar_file = covariate_file
  output:
    output_dir/"xqtl/summary_statistic"
  params:
    window_size = 50000
  resources:
    cpus_per_task = 8, mem = "32G", time = "8:00:00"
  shell:
    r"""
    set +eo

    mkdir -p $(dirname {params.out_prefix})
    apptainer exec -B /vol {singularity_image} qtltools \
      --vcf {input.vcf_file} \
      --bed {input.bed_file} \
      --out {output} \
    """
