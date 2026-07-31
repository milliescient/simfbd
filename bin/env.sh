# Shared toolchain paths (sourced by bin/run.sh). Edit here, not per study.
# Each is overridable from the environment, so a run can be pinned to another build
# without editing the repo. The manifest records the rb actually used.
RBIN=${RBIN:-/research/phyloworks/rb-sim/projects/cmake/build/rb}   # pinned sim build
RSCRIPT=${RSCRIPT:-/research/phyloworks/mm/root/envs/fsim/bin/Rscript}
PY=${PY:-/home/walker/.venv/bin/python}
