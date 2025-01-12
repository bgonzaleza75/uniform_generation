FROM "ghcr.io/rocker-org/devcontainer/r-ver:4.4" AS base

# Install additional packages
RUN install2.r --error --skipinstalled \
    ggplot2 \
    microbenchmark \
    tidyverse \
    rootSolve
