# syntax=docker/dockerfile:1.10
# RNAseq_STAR software in one image: STAR 2.7.11b (the official static binary, checksum as in
# 0_Setup/4_fetch_star.py), HTStream 1.3.3 (the Ceres module's exact build, in its own environment), SAMtools 1.17
# and Python 3.11 with requirements.txt, from the locked environments in pixi.lock. Base images are pinned by digest, so a rebuild of the same commit installs the same
# software. Built and published by .github/workflows/container.yml; on SCINet it runs with Apptainer
# (SOFTWARE_MODE=container).
FROM ghcr.io/prefix-dev/pixi:0.71.1-noble@sha256:2fe9f76c2baf6d68ebb939715828c1a02ceaf6033c14aa263c17d1ce5cf753bc AS build
WORKDIR /opt/rnaseq_star
COPY pixi.toml pixi.lock ./
RUN pixi install --locked --all
ADD --checksum=sha256:36e94b899a56b0ea5de5d65e722f55bc552b6713bd2d1e83a5c2abbd06c9881a \
    https://raw.githubusercontent.com/alexdobin/STAR/b1edc1208d91a53bf40ebae8669f71d50b994851/bin/Linux_x86_64_static/STAR \
    /opt/rnaseq_star/star/STAR
ADD --checksum=sha256:551e5bc4f91bf85b75c848aabba5acbcff4f4adf75e606cd17787038a155e69f \
    https://raw.githubusercontent.com/alexdobin/STAR/b1edc1208d91a53bf40ebae8669f71d50b994851/LICENSE \
    /opt/rnaseq_star/star/LICENSE
RUN chmod 0755 /opt/rnaseq_star/star/STAR

FROM ubuntu:24.04@sha256:24737cef17f4e22b1977df7d9a05e77c4e28b7da5b8911f8d8ccda136bf3573c
LABEL org.opencontainers.image.title="rnaseq_star" \
      org.opencontainers.image.description="RNAseq_STAR software: STAR 2.7.11b, HTStream 1.3.3, SAMtools 1.17, Python 3.11" \
      org.opencontainers.image.authors="Dawson Budke" \
      org.opencontainers.image.licenses="NOASSERTION"
COPY --from=build /opt/rnaseq_star/.pixi/envs/default /opt/rnaseq_star/.pixi/envs/default
COPY --from=build /opt/rnaseq_star/.pixi/envs/htstream /opt/rnaseq_star/.pixi/envs/htstream
COPY --from=build /opt/rnaseq_star/star /opt/rnaseq_star/star
ENV PATH=/opt/rnaseq_star/star:/opt/rnaseq_star/.pixi/envs/default/bin:/opt/rnaseq_star/.pixi/envs/htstream/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    CONDA_PREFIX=/opt/rnaseq_star/.pixi/envs/default \
    LC_ALL=C.UTF-8
RUN STAR --version && samtools --version | head -1 && hts_Stats --version && python3 -c "import pandas, numpy, scipy, Bio"
CMD ["bash"]
