# ADR #1: Implement VCF Chunking and Parallelization in Sash Workflow for PCGR Processing

**Status**: In Progress
**Date**: 2024-11-07
**Deciders**: Oliver Hofmann, Stephen Watts, Quentin Clayssen
**Technical Story**: Based on the limitations of PCGR in handling large variant datasets within the sash workflow, specifically impacting hypermutated samples.

## Context
[PCGR](https://sigven.github.io/pcgr/) (Personal Cancer Genome Reporter) currently has a variant processing limit of 500,000 variants per run. In the sash workflow, hypermutated samples often exceed this variant limit. PCGR has its own filtering steps, but an additional filtering step was also introduced in Bolt. By using VCF chunking and parallel processing, we can ensure that these large datasets are analyzed effectively without exceeding the PCGR variant limit, leading to larger annotation and a more scalable pipeline.

## Decision
To address the limitations of PCGR when handling hypermutated samples, we WILL implement the following:

1. **Split VCF Files into Chunks**: Input VCF files MUST be divided into chunks, each containing no more than 500,000 variants. This ensures that each chunk remains within PCGR’s processing capacity.

2. **Parallelize Processing**: Each chunk MUST be processed concurrently through PCGR to optimize processing time. The annotated outputs from all chunks MUST be merged to create a unified dataset.

3. **Integrate into Bolt Annotation**: The chunking and parallelization changes MUST be implemented in the Bolt annotation module to ensure seamless and scalable processing for large variant datasets.

4. **Efficiency Consideration**: For now, there MAY be a loss of efficiency for larger variant sets due to the fixed resources allocated for annotation. Further resource adjustments SHOULD be evaluated in the future.

## Consequences

### Positive Consequences
- **Improved Efficiency**: This approach allows large variant datasets to be processed within PCGR's constraints, enhancing efficiency and ensuring more comprehensive analysis.
- **Scalability**: Chunking and parallel processing make the sash workflow more scalable for hypermutated samples, accommodating larger datasets.

### Negative Consequences
- **Complexity**: Adding chunking and merging processes WILL increase complexity in data handling and ensuring integrity across all merged data.
- **Resource Demand**: Parallel processing MAY increase resource consumption, affecting system performance and requiring further resource management.

## Remaining Challenges
While the proposed approach mitigates the current limitations of PCGR, it MAY not fully resolve the issues for hypermutated samples with exceptionally high variant counts. Additional solutions MUST be explored, such as:

- **Additional Filtering Criteria**: Applying additional filters to reduce the variant count where applicable.
- **Alternative Reporting Methods**: Exploring more scalable reporting approaches that COULD handle higher variant loads.

## Status
**Status**: In Progress

## Links
- [Related PR for VCF Chunking and Parallelization Implementation](https://github.com/scwatts/bolt/pull/2)
- [PCGR Documentation on Variant Limit](https://sigven.github.io/pcgr/articles/running.html#large-input-sets-vcf)
- Discussion on Hypermutated Samples Handling
