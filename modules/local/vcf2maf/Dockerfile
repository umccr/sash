FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash coreutils procps findutils grep sed gawk \
    gzip bzip2 xz-utils tar zip unzip \
    curl wget ca-certificates git \
    perl libjson-perl libtext-csv-perl libarchive-zip-perl \
    libdbi-perl libdbd-sqlite3-perl \
    bcftools samtools tabix \
 && rm -rf /var/lib/apt/lists/*

# vcf2maf scripts
RUN git clone --depth 1 https://github.com/mskcc/vcf2maf.git /opt/vcf2maf \
 && ln -s /opt/vcf2maf/vcf2maf.pl /usr/local/bin/vcf2maf \
 && ln -s /opt/vcf2maf/maf2vcf.pl /usr/local/bin/maf2vcf