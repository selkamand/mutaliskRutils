test_that("mutalisk_to_dataframe works", {
  mutalisk_files <- dir(system.file("lusc_tcga", package = "mutaliskRutils"), pattern = "\\.txt$", full.names = TRUE)
  df_mutalisk_to_dataframe <- mutalisk_to_dataframe(mutalisk_files)
  df_expected <- structure(list(SampleID = c("TCGA-66-2789", "TCGA-66-2789", "TCGA-66-2789",
                                             "TCGA-66-2789", "TCGA-66-2789", "TCGA-66-2791", "TCGA-66-2791",
                                             "TCGA-66-2791", "TCGA-66-2791", "TCGA-66-2791", "TCGA-66-2792",
                                             "TCGA-66-2792", "TCGA-66-2792", "TCGA-66-2792", "TCGA-66-2793",
                                             "TCGA-66-2793", "TCGA-66-2794", "TCGA-66-2794", "TCGA-66-2794",
                                             "TCGA-66-2795", "TCGA-66-2795", "TCGA-66-2795", "TCGA-66-2795",
                                             "TCGA-66-2800", "TCGA-66-2800", "TCGA-66-2800", "TCGA-70-6722",
                                             "TCGA-70-6722", "TCGA-70-6722", "TCGA-70-6722", "TCGA-70-6723",
                                             "TCGA-70-6723", "TCGA-70-6723", "TCGA-85-6175", "TCGA-85-6175",
                                             "TCGA-85-6560", "TCGA-85-6560", "TCGA-85-6560", "TCGA-85-6560",
                                             "TCGA-85-6560", "TCGA-85-6561", "TCGA-85-6561", "TCGA-85-6561"
  ), Signatures = structure(c(1L, 3L, 4L, 5L, 2L, 1L, 3L, 4L, 5L,
                              2L, 1L, 3L, 4L, 2L, 4L, 6L, 3L, 4L, 6L, 1L, 3L, 4L, 2L, 3L, 4L,
                              6L, 1L, 3L, 4L, 2L, 4L, 5L, 6L, 1L, 5L, 1L, 3L, 4L, 6L, 2L, 3L,
                              4L, 2L), levels = c("1", "14", "2", "3", "5", "6"), class = "factor"),
  Contributions = c(0.05178, 0.0223, 0.7116, 0.16363, 0.05069,
                    0.05297, 0.07832, 0.68319, 0.09324, 0.09228, 0.13092, 0.02998,
                    0.71569, 0.12341, 0.89754, 0.10246, 0.0956, 0.69288, 0.21152,
                    0.05855, 0.02968, 0.7616, 0.15016, 0.04078, 0.75124, 0.20798,
                    0.0468, 0.05352, 0.81825, 0.08142, 0.70066, 0.13047, 0.16887,
                    0.26628, 0.73372, 0.06004, 0.01927, 0.7636, 0.0854, 0.07169,
                    0.02105, 0.82092, 0.15803)), row.names = c(NA, -43L), class = c("tbl_df",
                                                                                    "tbl", "data.frame"))

  expect_equal(df_mutalisk_to_dataframe, df_expected)
})

test_that("extract_sample_names_from_mutalisk_filenames works ",{
  mutalisk_files <- dir(system.file("lusc_tcga", package = "mutaliskRutils"), pattern = "\\.txt$", full.names = TRUE)
    expected <- c("TCGA-66-2789", "TCGA-66-2791", "TCGA-66-2792", "TCGA-66-2793",
      "TCGA-66-2794", "TCGA-66-2795", "TCGA-66-2800", "TCGA-70-6722",
      "TCGA-70-6723", "TCGA-85-6175", "TCGA-85-6560", "TCGA-85-6561")

    observed <- extract_sample_names_from_mutalisk_filenames(mutalisk_files)

    expect_equal(observed, expected)
})


test_that("extract_sample_names_from_mutalisk_files works",{
  mutalisk_files <- dir(system.file("lusc_tcga", package = "mutaliskRutils"), pattern = "\\.txt$", full.names = TRUE)
  expected <- c("TCGA-66-2789", "TCGA-66-2791", "TCGA-66-2792", "TCGA-66-2793",
                "TCGA-66-2794", "TCGA-66-2795", "TCGA-66-2800", "TCGA-70-6722",
                "TCGA-70-6723", "TCGA-85-6175", "TCGA-85-6560", "TCGA-85-6561")

  observed <- extract_sample_names_from_mutalisk_files(mutalisk_files)

  expect_equal(observed, expected)
})



test_that("mutalisk_to_dataframe works with sample_names from files", {
  mutalisk_files <- dir(system.file("lusc_tcga", package = "mutaliskRutils"), pattern = "\\.txt$", full.names = TRUE)
  df_mutalisk_to_dataframe_samplenames <- mutalisk_to_dataframe(mutalisk_files, sample_names_from_file_contents = TRUE)
  observed_sample_names <- sort(unique(df_mutalisk_to_dataframe_samplenames[['SampleID']]))

  sample_names <- extract_sample_names_from_mutalisk_filenames(mutalisk_files)
  expected_sample_names = sort(unique(sample_names))

  expect_equal(observed_sample_names, expected_sample_names)
})


test_that("mutalisk_best_signature_directory_to_dataframe works", {
  mutalisk_dir <- system.file("lusc_tcga", package = "mutaliskRutils")
  df_mutalisk_dir_to_df <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)
  df_expected <- structure(list(SampleID = c("TCGA-66-2789", "TCGA-66-2789", "TCGA-66-2789",
                                             "TCGA-66-2789", "TCGA-66-2789", "TCGA-66-2791", "TCGA-66-2791",
                                             "TCGA-66-2791", "TCGA-66-2791", "TCGA-66-2791", "TCGA-66-2792",
                                             "TCGA-66-2792", "TCGA-66-2792", "TCGA-66-2792", "TCGA-66-2793",
                                             "TCGA-66-2793", "TCGA-66-2794", "TCGA-66-2794", "TCGA-66-2794",
                                             "TCGA-66-2795", "TCGA-66-2795", "TCGA-66-2795", "TCGA-66-2795",
                                             "TCGA-66-2800", "TCGA-66-2800", "TCGA-66-2800", "TCGA-70-6722",
                                             "TCGA-70-6722", "TCGA-70-6722", "TCGA-70-6722", "TCGA-70-6723",
                                             "TCGA-70-6723", "TCGA-70-6723", "TCGA-85-6175", "TCGA-85-6175",
                                             "TCGA-85-6560", "TCGA-85-6560", "TCGA-85-6560", "TCGA-85-6560",
                                             "TCGA-85-6560", "TCGA-85-6561", "TCGA-85-6561", "TCGA-85-6561"
  ), Signatures = structure(c(1L, 3L, 4L, 5L, 2L, 1L, 3L, 4L, 5L,
                              2L, 1L, 3L, 4L, 2L, 4L, 6L, 3L, 4L, 6L, 1L, 3L, 4L, 2L, 3L, 4L,
                              6L, 1L, 3L, 4L, 2L, 4L, 5L, 6L, 1L, 5L, 1L, 3L, 4L, 6L, 2L, 3L,
                              4L, 2L), levels = c("1", "14", "2", "3", "5", "6"), class = "factor"),
  Contributions = c(0.05178, 0.0223, 0.7116, 0.16363, 0.05069,
                    0.05297, 0.07832, 0.68319, 0.09324, 0.09228, 0.13092, 0.02998,
                    0.71569, 0.12341, 0.89754, 0.10246, 0.0956, 0.69288, 0.21152,
                    0.05855, 0.02968, 0.7616, 0.15016, 0.04078, 0.75124, 0.20798,
                    0.0468, 0.05352, 0.81825, 0.08142, 0.70066, 0.13047, 0.16887,
                    0.26628, 0.73372, 0.06004, 0.01927, 0.7636, 0.0854, 0.07169,
                    0.02105, 0.82092, 0.15803)), row.names = c(NA, -43L), class = c("tbl_df",
                                                                                    "tbl", "data.frame"))
  expect_equal(df_mutalisk_dir_to_df, df_expected)
})

test_that("plot_stacked_bar works", {
  mutalisk_dir <- system.file("lusc_tcga", package = "mutaliskRutils")
  df_mutalisk <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)
  gg <- suppressMessages(plot_stacked_bar(df_mutalisk))
  vdiffr::expect_doppelganger(title = "Stacked Barplot", fig = gg)
})
