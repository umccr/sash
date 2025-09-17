FROM mambaorg/micromamba:1.5.8 AS build

USER root

RUN \
  micromamba create -y -p /env/ -c bioconda -c conda-forge \
    zlib && \
  micromamba clean --all --yes

FROM quay.io/biocontainers/hmftools-purple:4.2--hdfd78af_0

# NOTE(SW): copy over compatible build of BCFtools
COPY --from=quay.io/biocontainers/bcftools:1.19--h8b25389_0 /usr/local/ /usr/local/

COPY --from=build /env/ /env/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
