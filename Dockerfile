FROM --platform=linux/amd64 mambaorg/micromamba:1.0.0-bullseye-slim

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.yml /tmp/requirements.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /tmp/pyproject.toml
COPY --chown=$MAMBA_USER:$MAMBA_USER setup.cfg /tmp/setup.cfg
COPY --chown=$MAMBA_USER:$MAMBA_USER knn/ /tmp/knn/
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ /tmp/scripts/
RUN micromamba install -y -n base -c conda-forge -f /tmp/requirements.yml && \
        micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install /tmp/
