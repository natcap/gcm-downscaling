FROM --platform=linux/amd64 mambaorg/micromamba:1.0.0-bullseye-slim

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.yml /tmp/requirements.yml
RUN micromamba install -y -n base -c conda-forge -f /tmp/requirements.yml && \
        micromamba clean --all --yes
