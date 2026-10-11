# The tools the shell driver runs. simfbd sets all three before it calls the driver, from
# $RBIN and $RSCRIPT, else ~/.config/simfbd/tools.yml, else PATH; these are only what a driver
# called by hand falls back on.
RBIN=${RBIN:-rb}
RSCRIPT=${RSCRIPT:-Rscript}
PY=${PY:-python3}
