# First/last reporting: the extremes are always among the reported occurrences.
# Inherits every hyperprior from sims/batch1/config.sh.
MODEL=fbdr           # fbdr (dnFBDRMatrix) | bds | pyrate  -- only fbdr is wired
REPORTING=complete   # which specimens_<REPORTING>/ of the batch to read
INFER=complete       # dnFBDRMatrix complete=true
COND=time
RHO=1.0
GENS=200000; PRINTGEN=1000
NCORES=6
