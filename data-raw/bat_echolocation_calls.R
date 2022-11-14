## code to prepare `bat_echolocation_calls` dataset goes here
library(dplyr)

bat_echolocation_calls <-
  readRDS(
    "~/R/Projects/thesis/resevoir/harmonic_model/harmonic_model_dataset/traits.RDS"
  ) %>%
  mutate(
    call_recording = calls,
    harmonic_order = K,
    fundamental_frequency_contour = lapply(
      seq_along(f0_contour),
      function(i) cbind(t = time[[i]], f0 = f0_contour[[i]])
      )
  ) %>%
  select(
    bat,
    species,
    family,
    call_recording,
    fundamental_frequency_contour,
    duration,
    harmonic_order,
    dominant_harmonic
  )

usethis::use_data(bat_echolocation_calls, overwrite = TRUE)
