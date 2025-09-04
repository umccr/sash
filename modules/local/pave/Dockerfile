FROM --platform=linux/amd64 mambaorg/micromamba:1.5.7 as build

USER root

RUN \
  micromamba create -y -p /env/ -c bioconda -c conda-forge \
    zlib \
    bcftools=1.19 && \
  micromamba clean --all --yes

FROM quay.io/biocontainers/hmftools-pave:1.8--hdfd78af_0

COPY --from=build /env/ /env/

ENV PATH="/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/env/lib/:${LD_LIBRARY_PATH}"

# Ensure bcftools and companions are on default PATH expected by hmftools entrypoint
RUN ln -s /env/bin/bcftools /usr/local/bin/bcftools \
    && ln -s /env/bin/bgzip /usr/local/bin/bgzip \
    && ln -s /env/bin/tabix /usr/local/bin/tabix || true
