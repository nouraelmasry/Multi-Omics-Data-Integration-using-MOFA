library(MOFA2)
MOFAmodel <- load_model("model.hdf5")
plot_data_overview(MOFAmodel)
Nsamples = sum(MOFAmodel@dimensions$N)

sample_metadata <- data.frame(
  sample = samples_names(MOFAmodel)[[1]],
  condition = sample(c("A","B"), size = Nsamples, replace = T),
  age = sample(1:100, size = Nsamples, replace = T)
)

samples_metadata(MOFAmodel) <- sample_metadata
head(MOFAmodel@samples_metadata, n=3)

head(MOFAmodel@cache$variance_explained$r2_total[[1]])
head(MOFAmodel@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(MOFAmodel, x="view", y="factor")