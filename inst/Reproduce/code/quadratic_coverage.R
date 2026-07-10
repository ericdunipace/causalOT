#### Quadratic-truth coverage comparison for squared-IF variance ####

library(causalOT)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  as.integer(value)
}

env_int_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_num_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_chr_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  strsplit(value, ",", fixed = TRUE)[[1]]
}

ci_from_vcov <- function(effect, type, level = 0.95) {
  se <- sqrt(as.numeric(vcov(effect, type = type)))
  estimate <- as.numeric(coef(effect))
  crit <- stats::qnorm(1 - (1 - level) / 2)
  c(lwr = estimate - crit * se, upr = estimate + crit * se)
}

run_one <- function(seed, n, overlap, lambda, level) {
  set.seed(seed)
  data_gen <- Hainmueller$new(n = n, p = 6, design = "B", overlap = overlap)
  data_gen$gen_data()

  weights <- calc_weight(
    x = data_gen,
    estimand = "ATE",
    method = "COT",
    options = list(lambda = lambda)
  )
  effect <- estimate_effect(
    causalWeights = weights,
    model.function = stats::lm,
    augment.estimate = TRUE
  )

  estimate <- as.numeric(coef(effect))
  ci_if <- ci_from_vcov(effect, type = "if", level = level)
  ci_components <- ci_from_vcov(effect, type = "components", level = level)

  data.frame(
    seed = seed,
    n = n,
    overlap = overlap,
    lambda = lambda,
    estimate = estimate,
    variance_type = c("squared IF", "components"),
    ci_lwr = c(ci_if[["lwr"]], ci_components[["lwr"]]),
    ci_upr = c(ci_if[["upr"]], ci_components[["upr"]]),
    cover = c(
      ci_if[["lwr"]] <= 0 && ci_if[["upr"]] >= 0,
      ci_components[["lwr"]] <= 0 && ci_components[["upr"]] >= 0
    ),
    error = NA_character_
  )
}

nsims <- env_int("QUADRATIC_COVERAGE_NSIMS", 100L)
base_seed <- env_int("QUADRATIC_COVERAGE_SEED", 20260708L)
n_grid <- env_int_vec("QUADRATIC_COVERAGE_N", 512L)
overlap_grid <- env_chr_vec(
  "QUADRATIC_COVERAGE_OVERLAP",
  c("high", "medium", "low")
)
lambda_grid <- env_num_vec(
  "QUADRATIC_COVERAGE_LAMBDA",
  c(0, 0.01, 0.05, 0.1, 0.5, 1, 5)
)
level <- as.numeric(Sys.getenv("QUADRATIC_COVERAGE_LEVEL", unset = "0.95"))
output_file <- Sys.getenv(
  "QUADRATIC_COVERAGE_OUTPUT",
  unset = file.path("inst", "Reproduce", "data", "quadratic_coverage.rds")
)

grid <- expand.grid(
  sim = seq_len(nsims),
  n = n_grid,
  overlap = overlap_grid,
  lambda = lambda_grid,
  KEEP.OUT.ATTRS = FALSE
)
grid$seed <- base_seed + seq_len(nrow(grid))

results <- do.call(
  rbind,
  lapply(seq_len(nrow(grid)), function(i) {
    row <- grid[i, ]
    tryCatch(
      run_one(
        seed = row$seed,
        n = row$n,
        overlap = as.character(row$overlap),
        lambda = row$lambda,
        level = level
      ),
      error = function(e) {
        data.frame(
          seed = row$seed,
          n = row$n,
          overlap = as.character(row$overlap),
          lambda = row$lambda,
          estimate = NA_real_,
          variance_type = c("squared IF", "components"),
          ci_lwr = NA_real_,
          ci_upr = NA_real_,
          cover = NA,
          error = conditionMessage(e)
        )
      }
    )
  })
)

summary <- aggregate(
  cover ~ n + overlap + lambda + variance_type,
  data = results,
  FUN = function(x) mean(x, na.rm = TRUE)
)
names(summary)[names(summary) == "cover"] <- "coverage"
summary$n_success <- aggregate(
  cover ~ n + overlap + lambda + variance_type,
  data = results,
  FUN = function(x) sum(!is.na(x))
)$cover

note <- paste(
  "Appendix B.2 used a correctly specified linear outcome model, where the",
  "component and squared-IF variance estimators agree asymptotically, so those",
  "coverage figures should barely move. The meaningful finite-sample gap lives",
  "in the quadratic-truth settings behind Table 1. This experiment compares the",
  "legacy three-component variance with the squared-IF variance across the COT",
  "shrinkage grid; under-coverage by the component estimator in the shrinkage",
  "window is an empirical demonstration of the theory rather than an erratum."
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(results = results, summary = summary, note = note),
  file = output_file
)

print(summary)
cat(note, "\n")
