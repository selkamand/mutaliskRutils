test_that("mutalisk_to_dataframe works", {
  mutalisk_files <- dir(system.file("mutalisk_lusc_tcga",package = "mutalisk"), pattern = "\\.txt$", full.names = TRUE)
  df_mutalisk_to_dataframe <- mutalisk_to_dataframe(mutalisk_files)
  expect_snapshot(df_mutalisk)
})

test_that("mutalisk_to_dataframe works with sample_names specified", {
  mutalisk_files <- dir(system.file("mutalisk_lusc_tcga",package = "mutalisk"), pattern = "\\.txt$", full.names = TRUE)
  df_mutalisk_to_dataframe_samplenames <- mutalisk_to_dataframe(mutalisk_files)
  observed_sample_names <- sort(unique(df_mutalisk_to_dataframe_samplenames[['SampleID']]))

  sample_names <- mutalisk::extract_sample_names_from_mutalisk_filenames(mutalisk_files)
  expected_sample_names = sort(unique(sample_names))

  expect_equal(observed_sample_names, expected_sample_names)
})

test_that("mutalisk_best_signature_directory_to_dataframe works", {
  mutalisk_dir <- system.file("mutalisk_lusc_tcga",package = "mutalisk")
  df_mutalisk_dir_to_df <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)
  expect_snapshot(df_mutalisk_dir_to_df)
})

test_that("plot_stacked_bar works", {
  mutalisk_dir <- system.file("mutalisk_lusc_tcga",package = "mutalisk")
  df_mutalisk <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)
  gg <- suppressMessages(plot_stacked_bar(df_mutalisk))
  vdiffr::expect_doppelganger(title = "Stacked Barplot", fig = gg)
})
