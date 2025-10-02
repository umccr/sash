# umccr/sash: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

`sash` is the UMCCR post-processing WGS workflow. The workflow takes DRAGEN small variant calls and oncoanalyser results as input to perform annotation, prioritisation, rescue and filtering, and reporting for the WGS variant data. Additionally, `sash` runs several sensors for biomarker assessment and genomic characterisation including HRD status, mutational signatures, purity/ploidy, MSI, and TMB.

The general processes `sash` runs include:

- `gpgr` for generating the summary Cancer Report
- `PCGR` to report processed small somatic variants (annotated, rescued, filtered, prioritised)
- `CPSR` to report processed small germline variants (filtered, annotated, prioritised)
- `linxreport` to collate SV annotations and plots from LINX
- `MultiQC` for reporting various WGS statistics / metrics for QC
- `SAGE` variant calling to supplement DRAGEN small somatic variants
- `PURPLE` for TMB, MSI, CNV calling, and purity / ploidy estimation
- `HRDetect` and `CHORD` for HRD inference
- `MutationalPatterns` to fit mutational signatures
- `PAVE` for somatic variant annotation with MNV filtering

While the `sash` workflow utilises a range of tools and software, it is most closely coupled with [bolt](https://github.com/scwatts/bolt), a Python package that implements the UMCCR post-processing logic and supporting functionality.

- Requires Nextflow >= 22.10.6 and a container engine (Docker/Singularity/Apptainer/Podman).
- Uses GRCh38 reference data defined in `conf/refdata.config` accessed via `--ref_data_path`.
- Inputs are provided via a CSV samplesheet; no FASTQ inputs are expected by SASH.

## Requirements

- Java
- Nextflow ≥22.10.6
- A container engine, such as Docker, Singularity, Apptainer, or Podman

## Samplesheet input

Create a CSV samplesheet describing your DRAGEN and Oncoanalyser input directories. Pass it with `--input`. It must have a header row with columns: `id,subject_name,sample_name,filetype,filepath`.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

Provide one row per available input directory for a given analysis `id`. Rows sharing the same `id` must have the same `subject_name`. The `filetype` must be one of:

- `dragen_germline_dir`: directory containing DRAGEN germline outputs (normal sample)
- `dragen_somatic_dir`: directory containing DRAGEN somatic tumor/normal outputs (tumor sample)
- `oncoanalyser_dir`: directory containing HMFtools/Oncoanalyser outputs (e.g., Purple, Linx)

Example:

```csv
id,subject_name,sample_name,filetype,filepath
subject_a.example,subject_a,sample_germline,dragen_germline_dir,/path/to/dragen_germline/
subject_a.example,subject_a,sample_somatic,dragen_somatic_dir,/path/to/dragen_somatic/
subject_a.example,subject_a,sample_somatic,oncoanalyser_dir,/path/to/oncoanalyser/
```

Column descriptions:

| Column          | Description |
| --------------- | ----------- |
| `id`            | Analysis identifier used to group multiple rows belonging to the same subject/run. |
| `subject_name`  | Subject identifier. Must be identical across all rows for a given `id`. |
| `sample_name`   | Sample identifier. Used to set `normal_id` for `dragen_germline_dir` and `tumor_id` for `dragen_somatic_dir`. |
| `filetype`      | One of `dragen_germline_dir`, `dragen_somatic_dir`, `oncoanalyser_dir`. |
| `filepath`      | Absolute or relative path to the corresponding directory. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline and matches the format above.

### Required directory contents

SASH expects specific files inside each directory referenced in the samplesheet. Paths below are relative to each directory path you provide in `filepath`.

- `dragen_somatic_dir`
  - `${tumor_id}.hard-filtered.vcf.gz` and index `${tumor_id}.hard-filtered.vcf.gz.tbi`
  - `${tumor_id}.hrdscore.csv`
- `dragen_germline_dir`
  - `${normal_id}.hard-filtered.vcf.gz`
- `oncoanalyser_dir` (from nf-core/oncoanalyser)
  - `amber/` and `cobalt/` directories
  - `gridss/${tumor_id}.gridss.vcf.gz`
  - `sage/somatic/${tumor_id}.sage.somatic.vcf.gz` (+ `.tbi`)
  - `virusbreakend/` directory

Note: SASH runs PURPLE itself; precomputed PURPLE outputs are not required as inputs.

## Running the pipeline

Quickstart command:

```bash
nextflow run umccr/sash \
  --input samplesheet.csv \
  --ref_data_path /path/to/reference_data_root \
  --outdir results/ \
  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

This creates the following in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to reuse parameters across runs, specify them in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run umccr/sash -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: 'samplesheet.csv'
ref_data_path: '/path/to/reference_data_root'
outdir: 'results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Reference data

SASH reads all required resources from a single base directory passed via `--ref_data_path`. The internal paths and versions are defined in `conf/refdata.config` and include:

- Genome FASTA/FAI/DICT: `genomes/GRCh38_umccr/...`
- HMF reference data (WiGiTS): `hmf_reference_data/hmftools/<version>/...`
- PCGR bundle: `databases/pcgr/v<version>/`
- UMCCR resources: panels, known fusions, PoNs, config files

Organise your reference data root to match these relative paths, or provide a site config that overrides them. See `subworkflows/local/prepare_reference.nf` and `conf/refdata.config` for details.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull umccr/sash
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [umccr/sash releases page](https://github.com/umccr/sash/releases) and find the latest pipeline version (e.g. `0.5.0`). Then specify this when running with `-r` (single hyphen), for example `-r 0.5.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> 💡 If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

## Outputs

See detailed descriptions in `docs/output.md`. Top-level results include PCGR/CPSR HTML reports, cancer report, LINX/PURPLE artefacts, and a MultiQC summary.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
