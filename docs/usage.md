# umccr/sash: Usage

> Parameter documentation is generated automatically from `nextflow_schema.json`. Run `nextflow run umccr/sash --help`
> or use [nf-core/launch](https://nf-co.re/launch) for an interactive form.

## Introduction

umccr/sash is UMCCR’s post-processing pipeline for tumour/normal WGS analyses. It consumes DRAGEN secondary-analysis
outputs together with nf-core/oncoanalyser WiGiTS artefacts to perform small-variant rescue, annotation, filtering,
structural variant integration, PURPLE CNV calling, and reporting (PCGR, CPSR, GPGR cancer report, LINX, MultiQC).

- Requires Nextflow ≥ 22.10.6 and a container engine (Docker/Singularity/Apptainer/Podman/Conda).
- Uses GRCh38 reference data resolved from `--ref_data_path` (see [Reference data](#reference-data)).
- Expects inputs via a CSV samplesheet describing DRAGEN and Oncoanalyser directories; no FASTQ inputs are needed.

## Samplesheet input

Pass a CSV with `--input`. Each row represents one directory staged by upstream pipelines for a given analysis `id`.
Rows sharing the same `id` are grouped into a single tumour/normal run.

### Column definitions

| Column         | Description |
| -------------- | ----------- |
| `id`           | Unique analysis identifier grouping rows belonging to the same tumour/normal pair. |
| `subject_name` | Subject identifier; must be identical for all rows with the same `id`. |
| `sample_name`  | DRAGEN sample label. Used to derive tumour (`dragen_somatic_dir`) and normal (`dragen_germline_dir`) identifiers. |
| `filetype`     | One of the supported directory types below. |
| `filepath`     | Absolute or relative path to the directory containing the expected files. |

Example row set:

```csv
id,subject_name,sample_name,filetype,filepath
subject_a.example,subject_a,sample_germline,dragen_germline_dir,/path/to/dragen_germline/
subject_a.example,subject_a,sample_somatic,dragen_somatic_dir,/path/to/dragen_somatic/
subject_a.example,subject_a,sample_somatic,oncoanalyser_dir,/path/to/oncoanalyser/
```

An example sheet is included at `assets/samplesheet.csv`.

### Required directory contents

Paths below are relative to the value of `filepath` for each row. The pipeline targets nf-core/oncoanalyser ≥ 2.2.0
exports.

- `dragen_somatic_dir`
  - `<tumor_id>.hard-filtered.vcf.gz` and `<tumor_id>.hard-filtered.vcf.gz.tbi`
  - Optional: `<tumor_id>.hrdscore.csv` (ingested into the cancer report when present)
- `dragen_germline_dir`
  - `<normal_id>.hard-filtered.vcf.gz`
- `oncoanalyser_dir`
  - `amber/` and `cobalt/` directories (coverage inputs for PURPLE)
  - `sage_calling/somatic/<tumor_id>.sage.somatic.vcf.gz` (+ `.tbi`)
  - `esvee/<tumor_id>.esvee.ref_depth.vcf.gz` and accompanying directory (used to seed eSVee calling)
  - `virusbreakend/` directory

> SASH runs PURPLE internally; precomputed PURPLE outputs are not required as inputs.

## Running the pipeline

Quickstart command:

```bash
nextflow run umccr/sash \
  --input samplesheet.csv \
  --ref_data_path /path/to/reference_data_root \
  --outdir results/ \
  -profile docker
```

This launches the pipeline with the `docker` configuration profile. The following appear in the working directory:

```
work/           # Nextflow working files
results/        # Pipeline outputs (as specified by --outdir)
.nextflow_log   # Nextflow run log
```

### Parameter files and profiles

Reuse parameter sets via `-params-file params.yaml`:

```bash
nextflow run umccr/sash -profile docker -params-file params.yaml
```

```yaml
input: 'samplesheet.csv'
ref_data_path: '/data/refdata'
outdir: 'results/'
```

> ⚠️ Avoid using `-c` to pass pipeline parameters. `-c` should only point to Nextflow config files for resource tuning,
> executor settings or module overrides (see below).

You can generate YAML/JSON parameter files through [nf-core/launch](https://nf-co.re/launch) or Nextflow Tower.

## Reference data

All resources are resolved relative to `--ref_data_path` using `conf/refdata.config`. Confirm the directory contains the
expected subpaths (versions may change between releases):

- `genomes/GRCh38_umccr/` – GRCh38 FASTA, FAI and dict files plus sequence metadata.
- `hmf_reference_data/` – WiGiTS bundle with PURPLE GC profiles, eSVee panel-of-normals, SAGE hotspot resources, LINX
  transcripts and driver catalogues.
- `databases/pcgr/` – PCGR/CPSR annotation bundle.
- `umccr/` – bolt configuration files, driver panels, MultiQC templates, GPGR assets.
- `misc/` – panel-of-normals, APPRIS annotations, snpEff cache and other supporting data.

Refer to `docs/details.md` for a deeper breakdown of required artefacts.

## Nextflow configuration

### `-profile`

Profiles configure software packaging and cluster backends. Bundled profiles include `test`, `docker`, `singularity`,
`podman`, `shifter`, `charliecloud`, `apptainer` and `conda`. Combine multiple profiles with commas (later entries
override earlier ones). If no profile is supplied, Nextflow expects all software on `$PATH`, which is discouraged.

### `-resume`

Resume cached work by adding `-resume`. Nextflow matches stages using both file names and content; keep inputs identical
for cache hits. Supply a run name to resume a specific execution: `-resume <run-name>`. Use `nextflow log` to list
previous runs.

### `-c`

`-c custom.config` loads additional Nextflow configuration (eg. executor queues, resource overrides, institutional
profiles). See the [nf-core configuration docs](https://nf-co.re/docs/usage/configuration) for examples.

## Custom configuration

### Resource requests

Default resources suit typical datasets, but you can override CPUs/memory/time through custom config files. Many modules
honour nf-core’s automatic retry logic: certain exit codes trigger resubmission at 2× and 3× the original resources
before failing the run. Refer to the nf-core guides on
[max resources](https://nf-co.re/docs/usage/configuration#max-resources) and
[tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).

### Custom containers

nf-core pipelines default to Biocontainers/Bioconda images. You can override container or conda package selections in
config to use patched or institutional builds. See the
[updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section for patterns.

### Custom tool arguments

If you need to provide additional tool parameters beyond those exposed by pipeline options, set `process.ext.args`
(overrides per-module) or leverage module-specific hooks documented in nf-core. Review `conf/modules.config` for
supported overrides in umccr/sash.

## Outputs

See `docs/output.md` for a full description of generated artefacts (PCGR/CPSR HTML, cancer report, LINX, PURPLE, MultiQC
and supporting statistics).

## Running in the background

Nextflow supervises submitted jobs; keep the Nextflow process alive for the pipeline to finish. Options include:

- `nextflow run ... -bg` to launch detached and log to `.nextflow.log`.
- Using `screen`, `tmux` or similar to keep sessions alive.
- Submitting Nextflow itself through your scheduler (eg. `sbatch`), where it will launch child jobs.

## Nextflow memory requirements

The Nextflow JVM can request substantial RAM on large runs. Set an upper bound via environment variables, typically in
`~/.bashrc` or `~/.bash_profile`:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

Adjust limits to suit your environment.
