bootstrap: docker
from: mambaorg/micromamba:1.5.8-bookworm-slim

%files
  envs/environment.yml /opt

%post
  micromamba install -y -q -n base -f /opt/environment.yml
  micromamba remove -y -q -n base --prune quarto # quarto recipe has problems in this container https://github.com/conda-forge/quarto-feedstock/issues/30
  micromamba clean -ay
  apt-get update
  apt-get install procps wget -y
  wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.554/quarto-1.4.554-linux-amd64.deb -O /opt/quarto.deb
  apt-get install /opt/quarto.deb -y
  rm /opt/quarto.deb

%environment
  PATH=/opt/conda/bin:$PATH

%labels
  Author yang.e@wehi.edu.au
  Version v0.1
  XGBOOST_VERSION 1.6.2

%help
  This container is used in the MIBI pipeline to preprocess training data,
  train an XGBoost model, and apply the trained model.
