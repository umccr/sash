FROM mambaorg/micromamba:1.5.7 as build

USER root

RUN \
  micromamba create -y -p /env/ -c bioconda -c conda-forge \
    zlib && \
  micromamba clean --all --yes

FROM quay.io/biocontainers/hmftools-linx:2.1--hdfd78af_0

COPY --from=build /env/ /env/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"
