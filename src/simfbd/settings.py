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

    ("GAMMA_L_SHAPE", "rate_priors.gamma.speciation.shape",
     "PyRate's prior: an independent gamma draw per interval"),
    ("GAMMA_L_RATE", "rate_priors.gamma.speciation.rate",
     "left out, it is drawn per replicate from Uniform(0, estimated_rate_max) and estimated"),
    ("GAMMA_M_SHAPE", "rate_priors.gamma.extinction.shape", ""),
    ("GAMMA_M_RATE", "rate_priors.gamma.extinction.rate", ""),
    ("GAMMA_P_SHAPE", "rate_priors.gamma.fossil_sampling.shape", ""),
    ("GAMMA_P_RATE", "rate_priors.gamma.fossil_sampling.rate", "fixed; PyRate's is 1.1"),
    ("GAMMA_RATE_MAX", "rate_priors.gamma.estimated_rate_max",
     "upper bound of the uniform a speciation or extinction rate left out is drawn from"),

    ("RHO", "record.extant_sampling_probability",
     "chance a lineage alive at the present is seen there"),
    ("REPORT_PROB", "record.interior_reporting",
     "below 1, interior occurrences are thinned at a rate drawn per replicate from Uniform(0, 1)"),
    ("BIN_WIDTH", "record.occurrence_bins.width",
     "occurrence ages are reported to these bins"),
    ("BIN_MAX", "record.occurrence_bins.oldest_edge", ""),
    ("BIN_WIDTHS", "record.occurrence_bins.also_widths",
     "the same record on other grids too, as specimens_<model>_w<width>; 0 is the exact ages"),
    ("BIN_OFFSETS", "record.occurrence_bins.also_offsets",
     "each of those shifted by this fraction of a bin as well, as ..._o<offset>"),

    ("MIN_TAXA", "accepted_replicates.min_taxa",
     "a replicate outside these bounds is redrawn"),
    ("MAX_TAXA", "accepted_replicates.max_taxa", ""),
    ("TAXA_COUNT", "accepted_replicates.taxa_counted",
     "fossil: fossil-sampled taxa only; record: survivors seen at the present too"),
    ("MIN_FOSSIL_TAXA", "accepted_replicates.min_fossil_taxa", ""),
    ("N_EXTANT", "accepted_replicates.extant_species",
     "set, trees are drawn backward from exactly this many extant species and the origin is an "
     "output rather than a draw"),

    ("MAX_LINEAGES", "simulation_limits.max_lineages",
     "a draw expected to exceed this is refused unsimulated, which truncates the prior"),
    ("SIM_MAX_DRAWS", "simulation_limits.max_waiting_time_draws",
     "a draw whose simulation needs more waiting times than this is redrawn"),
]

RUN_SETTINGS = [
    ("MODEL", "inference.model", "fbdr: dnFBDRP; bds: dnBDS"),
    ("REPORTING", "inference.record_read", "the batch's specimens_<this> directory"),
    ("INFER", "inference.record_treated_as", "complete or incomplete"),
    ("AGES", "inference.occurrence_ages",
     "reported: the bins as they are; midpoint: each at its bin's middle; drawn: each uniform in its bin"),
    ("AGE_DRAWS", "inference.age_draws",
     "when drawn: this many draws, a chain each, pooled into one posterior per replicate"),
    ("COND", "inference.conditioned_on", ""),
    ("RHO", "inference.extant_sampling_probability", ""),
    ("ORIGIN_PRIOR", "inference.origin_prior", ""),
    ("SURVIVORS", "inference.survivors", ""),
    ("LATENTS", "inference.monitor_latents",
     "false drops the per-taxon b, d, tau_1 and tau_K from the log, which a rates-only "
     "comparison never reads and which costs 63 MB a replicate on a thousand ranges"),
    ("DROP_PRESENT_ONLY", "inference.drop_present_only_taxa", ""),
    ("GENS", "mcmc.generations", ""),
    ("PRINTGEN", "mcmc.sample_every", ""),
    ("SEED", "mcmc.seed", "each replicate's chain is seeded with this plus its number"),
    ("NCORES", "mcmc.chains_at_once", "replicates analyzed in parallel"),
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


# what a batch cannot be simulated without, and a run cannot be analyzed without
BATCH_REQUIRED = ["NINTERVALS", "INTERVAL_WIDTH", "GAMMA_L_SHAPE", "GAMMA_M_SHAPE", "GAMMA_P_SHAPE",
                  "GAMMA_P_RATE", "AGE_MIN", "AGE_MAX", "BIN_WIDTH", "BIN_MAX", "ORIGIN_SAMPLED"]
RUN_REQUIRED = ["MODEL", "REPORTING", "INFER", "COND", "RHO", "GENS", "PRINTGEN", "NCORES"]
# The value a setting takes when a config leaves it out is in defaults/batch.yml and
# defaults/run.yml, not here.

# a setting that takes one of a few words
CHOICES = {"MODEL": ("fbdr", "bds"), "INFER": ("complete", "incomplete"),
           "TAXA_COUNT": ("fossil", "record"),
           "ORIGIN_PRIOR": ("uniform", "exponential", "none"),
           "AGES": ("reported", "midpoint", "drawn")}
