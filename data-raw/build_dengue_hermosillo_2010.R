source("R/model_utils.R", local = FALSE)

raw_data <- load_study_data(
  path = package_path(
    "extdata",
    "raw",
    "study_data_dengue_hemorrhagic_fever.csv"
  )
)

dengue_hermosillo_2010 <- raw_data

save(
  dengue_hermosillo_2010,
  file = "data/dengue_hermosillo_2010.rda",
  compress = "bzip2"
)
