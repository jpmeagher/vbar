cavi_plvm <- function(
  plvm_list,
  tol = 1e-6, max_iter = 1000,
  n_samples = 1000, random_seed = NULL,
  progress_bar = TRUE,
  perform_checks = TRUE
){
  elbo <- c(
    -Inf,
    compute_plvm_elbo(
      plvm_list,
      n_samples = n_samples, random_seed = random_seed,
      perform_checks = perform_checks
    )
  )
  i <- 1
  if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = max_iter, style = 3)
  while (i <= max_iter) {
    elbo <- c(
      elbo,
      compute_plvm_elbo(
        plvm_list,
        n_samples = n_samples, random_seed = random_seed,
        perform_checks = FALSE
      )
    )
    if (progress_bar) utils::setTxtProgressBar(pb, i)
    i <- i + 1
  }
  if (progress_bar) close(pb)

  elbo
}
