# Shared toolchain paths (sourced by bin/run.sh). Edit here, not per study.
# Each is overridable from the environment, so a run can be pinned to another build
# without editing the repo. The manifest records the rb actually used.
# fbdr merged with paleo-rate-estimators: the one binary every arm of the final SBC runs
# on, because no earlier build has both dnBDS(origin_prior=) and the estimator functions.
# It lives in a worktree, which must not be removed while the study is running.
RBIN=${RBIN:-/research/phyloworks/rb-fbdr-est/projects/cmake/build-fbdr-est/rb}
RSCRIPT=${RSCRIPT:-/research/phyloworks/mm/root/envs/fsim/bin/Rscript}
PY=${PY:-/home/walker/.venv/bin/python}
