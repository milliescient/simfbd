"""Every setting a batch or a run has: the name it is written under, and what it means.

A config is what you ask for, in these names, and a manifest is the record of what was done. The first column is the name the simulator and the Rev templates know a setting by,
which is how it reaches them and nothing a person writes or reads.
"""

# internal name, the name it is written under, what it means
BATCH_SETTINGS = [
    ("NINTERVALS", "timeline.rate_intervals",
     "piecewise-constant rate intervals"),
    ("INTERVAL_WIDTH", "timeline.interval_width",
     "width of each interval; the oldest runs on to the origin"),
    ("PRESENT", "timeline.record_stops_at",
     "the record ends here, not at zero; younger occurrences are discarded"),
    ("AGE_MIN", "origin_age.min", "origin ~ Uniform(min, max)"),
    ("AGE_MAX", "origin_age.max", ""),
    ("ORIGIN_SAMPLED", "origin_age.origin_lineage_must_be_sampled", ""),

    ("LMEAN", "rate_priors.speciation.log_mean",
     "lognormal, on the log scale; one draw per interval"),
    ("LSD", "rate_priors.speciation.log_sd", ""),
    ("MMEAN", "rate_priors.extinction.log_mean", ""),
    ("MSD", "rate_priors.extinction.log_sd", ""),
    ("TMEAN", "rate_priors.turnover.log_mean",
     "stationary prior only: extinction is speciation times turnover"),
    ("TSD", "rate_priors.turnover.log_sd", ""),
    ("PMEAN", "rate_priors.fossil_sampling.log_mean", ""),
    ("PSD", "rate_priors.fossil_sampling.log_sd", ""),
    ("LAMBDA_A", "rate_priors.anagenetic_speciation",
     "above 0, an anagenetic rate is drawn from the speciation prior; 0 is budding only"),

    ("GMRF_SD", "rates_through_time.random_walk.step_sd",
     "sd of the step between adjacent intervals' log rates"),
    ("GMRF_SD_MEAN", "rates_through_time.random_walk.step_sd_prior.log_mean",
     "speciation and extinction: the step sd is drawn per replicate, lognormal"),
    ("GMRF_SD_SD", "rates_through_time.random_walk.step_sd_prior.log_sd", ""),
    ("GMRF_PSI_MEAN", "rates_through_time.random_walk.fossil_sampling_step_sd_prior.log_mean",
     "fossil sampling gets its own step sd"),
    ("GMRF_PSI_SD", "rates_through_time.random_walk.fossil_sampling_step_sd_prior.log_sd", ""),
    ("OU_THETA", "rates_through_time.stationary.mean_reversion",
     "fixed; set this or correlation_prior, not both"),
    ("OU_DREF", "rates_through_time.stationary.reference_lag",
     "the lag the correlation is stated at"),
    ("OU_PHI_A", "rates_through_time.stationary.correlation_prior.beta_a",
     "correlation at the reference lag, drawn per replicate from Beta(a, b)"),
    ("OU_PHI_B", "rates_through_time.stationary.correlation_prior.beta_b", ""),
    ("LANCH", "rates_through_time.stationary.level_sd.speciation",
     "sd of each replicate's mean log rate about the prior mean; 0 fixes it there"),
    ("TANCH", "rates_through_time.stationary.level_sd.turnover", ""),
    ("PANCH", "rates_through_time.stationary.level_sd.fossil_sampling", ""),
    ("OU_SD_MEAN", "rates_through_time.stationary.spread_prior.log_mean",
     "speciation and turnover: the sd through time is drawn per replicate, lognormal"),
    ("OU_SD_SD", "rates_through_time.stationary.spread_prior.log_sd", ""),
    ("OU_PSI_MEAN", "rates_through_time.stationary.fossil_sampling_spread_prior.log_mean",
     "fossil sampling gets its own"),
    ("OU_PSI_SD", "rates_through_time.stationary.fossil_sampling_spread_prior.log_sd", ""),

    ("RHO", "record.extant_sampling_probability",
     "chance a lineage alive at the present is seen there"),
    ("REPORT_PROB", "record.interior_reporting",
     "below 1, interior occurrences are thinned at a rate drawn per replicate from Uniform(0, 1)"),
    ("BIN_WIDTH", "record.occurrence_bins.width",
     "occurrence ages are reported to these bins"),
    ("BIN_MAX", "record.occurrence_bins.oldest_edge", ""),

    ("MIN_TAXA", "accepted_replicates.min_taxa",
     "a replicate outside these bounds is redrawn"),
    ("MAX_TAXA", "accepted_replicates.max_taxa", ""),
    ("TAXA_COUNT", "accepted_replicates.taxa_counted",
     "fossil: fossil-sampled taxa only; record: survivors seen at the present too"),
    ("MIN_FOSSIL_TAXA", "accepted_replicates.min_fossil_taxa", ""),

    ("MAX_LINEAGES", "simulation_limits.max_lineages",
     "a draw expected to exceed this is refused unsimulated, which truncates the prior"),
    ("SIM_MAX_DRAWS", "simulation_limits.max_waiting_time_draws",
     "a draw whose simulation needs more waiting times than this is redrawn"),
]

RUN_SETTINGS = [
    ("MODEL", "inference.model", "fbdr: dnFBDRP; bds: dnBDS; fbdsp: dnFBDSP"),
    ("REPORTING", "inference.record_read", "the batch's specimens_<this> directory"),
    ("INFER", "inference.record_treated_as", "complete or incomplete"),
    ("COND", "inference.conditioned_on", ""),
    ("RHO", "inference.extant_sampling_probability", ""),
    ("ORIGIN_PRIOR", "inference.origin_prior", ""),
    ("SURVIVORS", "inference.survivors", ""),
    ("DROP_PRESENT_ONLY", "inference.drop_present_only_taxa", ""),
    ("GENS", "mcmc.generations", ""),
    ("PRINTGEN", "mcmc.sample_every", ""),
    ("NCORES", "mcmc.chains_at_once", "replicates fitted in parallel"),
]

# what a manifest records beyond settings: the name a script asks for, and where it sits
BATCH_RECORD = {
    "kind": "kind", "batch": "batch", "config_hash": "config_hash",
    "true_vals_order": "true_vals_order", "NREPS": "replicates",
    "filtered": "prior_filter.refused", "filtered_of": "prior_filter.of",
    "filtered_frac": "prior_filter.fraction",
}
RUN_RECORD = {
    "kind": "kind", "run": "run", "batch": "batch",
    "batch_config_hash": "batch_config_hash", "config_hash": "config_hash",
    "true_vals_order": "true_vals_order", "NREPS": "replicates", "specimens": "specimens",
    "infer_script_hash": "toolchain.infer_script_hash",
    "rb_path": "toolchain.rb_path", "rb_md5": "toolchain.rb_md5",
    "rb_describe": "toolchain.rb_describe",
    # the names a manifest.tsv used for these two
    "survivors": "inference.survivors", "drop_present_only": "inference.drop_present_only_taxa",
}


# what a batch cannot be simulated without, and a run cannot be fitted without
BATCH_REQUIRED = ["NINTERVALS", "INTERVAL_WIDTH", "LMEAN", "LSD", "MMEAN", "MSD", "PMEAN", "PSD",
                  "AGE_MIN", "AGE_MAX", "BIN_WIDTH", "BIN_MAX", "LAMBDA_A", "ORIGIN_SAMPLED"]
RUN_REQUIRED = ["MODEL", "REPORTING", "INFER", "COND", "RHO", "GENS", "PRINTGEN", "NCORES"]
# The value a setting takes when a config leaves it out is in defaults/batch.yml and
# defaults/run.yml, not here.

# a setting that takes one of a few words
CHOICES = {"MODEL": ("fbdr", "bds", "fbdsp"), "INFER": ("complete", "incomplete"),
           "TAXA_COUNT": ("fossil", "record")}
