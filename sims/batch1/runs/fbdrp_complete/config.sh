# Complete reporting: every fossil kept at its exact age.
# Inherits every hyperprior from sims/batch1/config.sh.
MODEL=fbdr           # fbdr (dnFBDRMatrix) | bds | pyrate  -- only fbdr is wired
REPORTING=complete   # which specimens_<REPORTING>/ of the batch to read
INFER=complete       # dnFBDRMatrix complete=true (every fossil kept)
COND=time
RHO=1.0
GENS=10000000; PRINTGEN=10000
NCORES=16
