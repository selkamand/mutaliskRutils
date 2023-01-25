test_that("mutalisk_to_dataframe works", {
  mutalisk_files <- dir(system.file("lusc_tcga",package = "mutalisk"), pattern = "\\.txt$", full.names = TRUE)
  df_mutalisk_to_dataframe <- mutalisk_to_dataframe(mutalisk_files)
  expect_snapshot(df_mutalisk_to_dataframe)
})

test_that("extract_sample_names_from_mutalisk_filenames works ",{
  mutalisk_files <- dir(system.file("lusc_tcga",package = "mutalisk"), pattern = "\\.txt$", full.names = TRUE)
    expected <- c("TCGA-66-2789", "TCGA-66-2791", "TCGA-66-2792", "TCGA-66-2793",
      "TCGA-66-2794", "TCGA-66-2795", "TCGA-66-2800", "TCGA-70-6722",
      "TCGA-70-6723", "TCGA-85-6175", "TCGA-85-6560", "TCGA-85-6561")

    observed <- extract_sample_names_from_mutalisk_filenames(mutalisk_files)

    expect_equal(observed, expected)
})


test_that("extract_sample_names_from_mutalisk_files works",{
  mutalisk_files <- dir(system.file("lusc_tcga",package = "mutalisk"), pattern = "\\.txt$", full.names = TRUE)
  expected <- c("TCGA-66-2789", "TCGA-66-2791", "TCGA-66-2792", "TCGA-66-2793",
                "TCGA-66-2794", "TCGA-66-2795", "TCGA-66-2800", "TCGA-70-6722",
                "TCGA-70-6723", "TCGA-85-6175", "TCGA-85-6560", "TCGA-85-6561")

  observed <- extract_sample_names_from_mutalisk_files(mutalisk_files)

  expect_equal(observed, expected)
})



test_that("mutalisk_to_dataframe works with sample_names from files", {
  mutalisk_files <- dir(system.file("lusc_tcga",package = "mutalisk"), pattern = "\\.txt$", full.names = TRUE)
  df_mutalisk_to_dataframe_samplenames <- mutalisk_to_dataframe(mutalisk_files, sample_names_from_file_contents = TRUE)
  observed_sample_names <- sort(unique(df_mutalisk_to_dataframe_samplenames[['SampleID']]))

  sample_names <- extract_sample_names_from_mutalisk_filenames(mutalisk_files)
  expected_sample_names = sort(unique(sample_names))

  expect_equal(observed_sample_names, expected_sample_names)
})


test_that("mutalisk_best_signature_directory_to_dataframe works", {
  mutalisk_dir <- system.file("lusc_tcga",package = "mutalisk")
  df_mutalisk_dir_to_df <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)
  expect_snapshot(df_mutalisk_dir_to_df)
})

test_that("plot_stacked_bar works", {
  mutalisk_dir <- system.file("lusc_tcga",package = "mutalisk")
  df_mutalisk <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)
  gg <- suppressMessages(plot_stacked_bar(df_mutalisk))
  vdiffr::expect_doppelganger(title = "Stacked Barplot", fig = gg)
})
