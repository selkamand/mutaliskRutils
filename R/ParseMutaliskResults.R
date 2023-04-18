# declarations
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Mutalisk to dataframe
#'
#' You probably want to use tags$strong{mutalisk_to_dataframe} instead.
#' See ?mutalisk_to_dataframe
#'
#' @param mutalisk_file a vector of filepaths, each leading to the report.txt files output when downloading best_signature results for all vcfs in cohort
#' @param sample_names_from_file_contents  guess sample names from filecontents instead of filenames (flag)
#'
#' @return tibble
#' @importFrom rlang .data
mutalisk_to_dataframe_single_sample <- function(mutalisk_file, sample_names_from_file_contents = FALSE){
  mutfile_v = readLines(mutalisk_file)

  #browser()
  #Ensure file is mutalisk file
  assertions::assert_identical(mutfile_v[1],y = "#####SIGNATURE DECOMPOSITION REPORT", msg = "Invalid mutalisk file. First line should be `#####SIGNATURE DECOMPOSITION REPORT`, not {arg_value}")

  sig_numbers_s = mutfile_v[grep(pattern = "\\*\\*THIS_SIGs ", x = mutfile_v)][1]
  sig_numbers_s = sig_numbers_s %>% strsplit(" ") %>% unlist()
  sig_numbers_f = sig_numbers_s[-1] %>% as.factor()

  sig_contributions_s = mutfile_v[grep(pattern = "\\*\\*THIS_SIG_CONTRIBUTIONS ", x = mutfile_v)][1]
  sig_contributions_s = sig_contributions_s %>% strsplit(" ") %>% unlist()
  sig_contributions_n = sig_contributions_s[-1] %>% as.numeric()

  if(!sample_names_from_file_contents)
    patient_id_s = extract_sample_names_from_mutalisk_filenames(mutalisk_file)
  else
    patient_id_s = extract_sample_names_from_mutalisk_files(mutalisk_file)


  data.frame(SampleID = rep(patient_id_s, times = length(sig_numbers_f)), Signatures = sig_numbers_f, Contributions =  sig_contributions_n) %>%
    dplyr::tibble()
}

#' Mutalisk files to dataframe
#'
#' @param mutalisk_files a vector of filepaths, each leading to the report.txt files output when downloading best_signature results for all vcfs in cohort
#' @param sample_names_from_file_contents guess sample names from filecontents instead of filenames (flag)
#' @return  a dataframe containing three columns:
#' \enumerate{
#' \item \strong{SampleID}: a sample identifier.
#' \item \strong{Signatures}: an identifier for a particular signature.
#' \item \strong{Contributions}: the percentage contribution of the signature to the patients genetic profile (0.1 = 10 percent).
#' }
#' @export
#'
mutalisk_to_dataframe <- function(mutalisk_files, sample_names_from_file_contents = FALSE){
  assertions::assert_file_exists(mutalisk_files)
  assertions::assert_flag(sample_names_from_file_contents)

  mutalisk_ls <- purrr::map(
    seq_along(mutalisk_files),
    function(i){
      mutalisk_to_dataframe_single_sample(
        mutalisk_file = mutalisk_files[i],
        sample_names_from_file_contents = sample_names_from_file_contents
        )
    }
  )
  mutalisk_df <- purrr::list_rbind(mutalisk_ls)


  return(mutalisk_df)
}

#' Sample Names From Mutalisk Output
#' @param mutalisk_filenames names of mutalisk output files (character)
#'
#' @return sample name (string)
extract_sample_names_from_mutalisk_filenames  <- function(mutalisk_filenames){
  assertions::assert_character(mutalisk_filenames)

  basename(mutalisk_filenames) %>%
    sub(pattern = ".*mutalisk_input_(.*?)\\..*$", replacement = "\\1", x =  mutalisk_filenames) %>%
    sub(pattern = "(.*)_.*$", replacement = "\\1", x = .) %>%
    return()
}

#' Sample Names From Mutalisk File Contents
#' @param mutalisk_filenames names of mutalisk output files (character)
#'
#' @return sample name (string)
extract_sample_names_from_mutalisk_files <- function(mutalisk_filenames){
  ###./uploads/mutalisk_input_TCGA-66-2789.HcQhu3ih8e_whole_sig.txt

  sample_names <- vapply(
    mutalisk_filenames,

    function(filename){
      secondline <- readLines(filename, n = 2)[[2]]
      sample_name <- secondline %>%
        sub(pattern = ".*mutalisk_input_(.*?)\\/*$", replacement = "\\1", x =  .) %>%
        sub(pattern = "(.*?)\\..*$", replacement = "\\1", x = .)
      return(sample_name)
  }, FUN.VALUE = character(1))

  unname(sample_names)
}

#' Mutalisk directory to dataframe
#'
#'
#' @param directory path to \strong{mutalisk_best_fit} folder.
#' To obtain, run your VCFs through mutalisk. Select Mutational Signature (Best only) and click 'Get the selected result for all samples at once'.
#' Then unzip the file, and youre ready to go
#' @inheritParams mutalisk_dataframe_inform_user_of_metadata
#'
#' @return tibble
#' @export
#'
mutalisk_best_signature_directory_to_dataframe <- function(directory, metadata = NA){
 assertions::assert_directory_exists(directory)

  filenames = dir(full.names = TRUE, path = directory, pattern = "mutalisk_input.*\\.txt$")

  df = mutalisk_to_dataframe(filenames)

  if(!is.na(metadata)){
    df = mutalisk_dataframe_inform_user_of_metadata(df, metadata)
  }

  if(nrow(df) == 0) {
    stop("Failed to find mutalisk files inside directory. Are you sure this is the appropriate mutalisk output?")
  }


  return(df)
}


# Plots -------------------------------------------------------------------

#' Cohort-Level Mutational Signature Visualisation
#'
#' A note of warning: for different mutalisk runs, this function will not enforce uniform colours for a single mutational signature. Better to get all data in at once and add a facet_wrap call
#' @importFrom pals kovesi.diverging_rainbow_bgymr_45_85_c67
#'
#' @param mutalisk_dataframe a dataframe that can be produced using mutalisk_best_signature_directory_to_dataframe. Can also just make it yourself, if you want to visualise non-mutalisk data. Dataframe just needs 3 columns:
#' \enumerate{
#' \item \strong{SampleID}: a sample identifier.
#' \item \strong{Signatures}: an identifier for a particular signature.
#' \item \strong{Contributions}: the percentage contribution of the signature to the patients genetic profile (0.1 = 10%).
#' }
#' @param lump_type one of "min_prop", "topn", or "none".
#'
#' \strong{min_prop} will allow lump together all signatures that contribute less than \strong{lump_min}.
#' \strong{topn} will keep the topn contributing signatures distinct, and lump the rest together (string)
#' @param lump_min see \strong{lump_type}
#' @param topn see \strong{lump_type}
#' @param legend Where should the legend be placed? ("top", "left", "bottom", "right")
#' @param legend_direction How should the legend be oriented. By defualt will guess based on position of legend ("vertical", "horizontal")
#' @param pal Palette to use for generating colours.
#' @param facet_column name of column to use for faceting (string)
#' @param fontsize_strip fontsize of facet titles (number)
#' @param color_of_other colour of signatures lumped into 'other' (string)
#' @param fontsize_axis_title fontsize of axis titles (number)
#'
#' @return a ggplot (gg)
#' @export
plot_stacked_bar <- function(mutalisk_dataframe, lump_type = "min_prop", lump_min=0.1,
                             topn = 5, legend="right",
                             legend_direction = NA, pal = pals::kovesi.diverging_rainbow_bgymr_45_85_c67,
                             color_of_other = "grey60", facet_column = NA, fontsize_strip = 18,
                             fontsize_axis_title = 18){
  assertions::assert_names_include(mutalisk_dataframe, names = c("Signatures", "SampleID", "Contributions"))
  checkmate::assert_choice(x = lump_type, choices = c("min_prop", "topn", "none"))
  checkmate::assert_choice(x = legend, choices = c("top", "left", "bottom", "right"))
  assertions::assert_number(fontsize_strip)
  assertions::assert_number(fontsize_axis_title)

  if(lump_type != "min_prop") stop("Other lump_types are not currently supported")

  if(is.na(legend_direction))
    legend_direction = ifelse(legend=="right" | legend == "left", yes = "vertical", no = "horizontal")
  else
    checkmate::assert_choice(legend_direction, choices = c("horizontal", "vertical"))

  #Ungroup Variables
  mutalisk_dataframe <- mutalisk_dataframe %>%
    dplyr::ungroup()

  if (lump_type == "min_prop") {
    assertions::assert_number(x = lump_min)

    message("Lumping together signatures with contributions < ", lump_min, " in all samples as 'Other'")
    mutalisk_dataframe <- mutalisk_dataframe %>%
      dplyr::group_by(.data[["Signatures"]]) %>%
      dplyr::mutate(SignaturesAllContributionsBelowMin = all(.data[["Contributions"]] < lump_min)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Signatures = ifelse(.data[['SignaturesAllContributionsBelowMin']], yes="Other", no = as.character(.data[['Signatures']])))

    #return(mutalisk_dataframe)
    mutalisk_dataframe$Signatures <- mutalisk_dataframe$Signatures %>%
      forcats::fct_relevel("Other", after = Inf) %>%
      droplevels()
  }
  else(
    stop("lump type: ", lump_type, "is not yet supported")
    )

  # if (lump_type == "topn") {
  #   checkmate::assert_number(x = topn, finite = TRUE)
  #   message(
  #     "Lumping together all signatures except those with the ",
  #     topn ,
  #     " greatest contributions"
  #   )
  #   mutalisk_dataframe <- mutalisk_dataframe %>%
  #     dplyr::group_by(SampleID) %>%
  #     dplyr::mutate(
  #       Signatures = forcats::fct_lump_n(Signatures, w = Contributions, n = topn)
  #     )
  # }
  #
  # if (lump_type == "none") {
  #   message("No signatures will be lumped together as 'other'. All signatures in a mutalisk 'best fit' set will be shown.")
  # }

  #Reorder signature levels by contribution size:
  mutalisk_dataframe <- mutalisk_dataframe %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Signatures = forcats::fct_reorder(.data[['Signatures']], .data[['Contributions']], .desc = TRUE))

  if ("Other" %in% levels(mutalisk_dataframe$Signatures)){
    mutalisk_dataframe$Signatures <- mutalisk_dataframe$Signatures %>%
      forcats::fct_relevel("Other", after = Inf) %>%
      droplevels()
  }

  gg = mutalisk_dataframe %>% ggplot2::ggplot() +
    ggplot2::geom_col(
      ggplot2::aes(
        y = .data[['SampleID']],
        x = .data[['Contributions']],
        fill = .data[['Signatures']],
      ), position = "fill") +
    ggplot2::ylab("Sample") +
    ggplot2::xlab("Contribution") +
    ggplot2::scale_fill_manual(
      values = c(
        pal(n = nlevels(mutalisk_dataframe$Signatures) - 1),color_of_other
      )
    ) +
    ggplot2::labs(fill = "Signature") +
    ggthemes::theme_fivethirtyeight() +
    ggplot2::theme(legend.position = legend, legend.direction = "vertical") +
    ggplot2::theme(axis.title = ggplot2::element_text(face = "bold", size = fontsize_axis_title)) +
    ggplot2::theme(strip.text= ggplot2::element_text(face = "bold", size = fontsize_strip), strip.background = ggplot2::element_blank()) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey80")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(), plot.background = ggplot2::element_blank())

  if(!is.na(facet_column)){
    checkmate::assert_choice(facet_column, choices = names(mutalisk_dataframe))
    gg <- gg + ggplot2::facet_grid(rows  = ggplot2::vars(!!dplyr::sym(facet_column)), scales = "free", space = "free")
  }

  mutalisk_dataframe_metadata_column_message(mutalisk_dataframe)

  return(gg)
}

#' Plot Signature-Level Dotplot
#' Plots a signature-Level dotplot
#'
#' @inheritParams plot_stacked_bar
#'
#' @return a ggplot object
#' @export
#'
plot_signature_contribution_jitterplot <- function(mutalisk_dataframe){
  requireNamespace("stats", quietly = TRUE)

  assertions::assert_names_include(mutalisk_dataframe, c("Signatures", "SampleID", "Contributions"))

  levels(mutalisk_dataframe$Signatures) <- gtools::mixedsort(levels((mutalisk_dataframe$Signatures)))

  mutalisk_dataframe <- mutalisk_dataframe_expand(mutalisk_dataframe)

  mutalisk_dataframe_metadata_column_message(mutalisk_dataframe)

  mutalisk_dataframe %>%
    ggplot2::ggplot(ggplot2::aes(.data[['Signatures']], .data[['Contributions']])) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::stat_summary(
      fun = stats::median, fun.min = stats::median, fun.max = stats::median,
      geom = "crossbar", color = "red", width = 0.7, lwd = 0.4) +
    ggplot2::theme_bw()
}

#' Expand mutalisk_dataframe
#'
#' In normal mutalisk dataframe, each sample has data for ONLY the 1-7 signatures that comprise the 'best fit'.
#' Bascically, we end up with a dataframe where not all signatures have entries for all samples. We call this 'implicit' missing values.
#' This means that a signature level jitterplot won't show the samples where it was not included in this 'best fit' set, when we'd actually want to know that it contributed 0% to that sample.
#' This function fixes the issue by adding entries for ALL signature - sample pairs, with Contributions set to 0% where relevant.
#'
#' @inheritParams plot_stacked_bar
#'
#' @return dataframe containing all combinations of Sample ID and Signatures.
#' For cases where a signature was not included in the 'best fit' subset, Contribution is set to 0%.
#'
mutalisk_dataframe_expand <- function(mutalisk_dataframe){
  assertions::assert_names_include(mutalisk_dataframe, names = c("Signatures", "SampleID", "Contributions"))

  #Make implicit missing values explicit
  mutalisk_dataframe %>%
    tidyr::complete(.data[['Signatures']], .data[['SampleID']], fill = list(Contributions = 0)) %>%
    return()
}


# Metadata ----------------------------------------------------------------
#' Add sample metadata
#'
#' Add sample metadata from dataframe to mutalisk data.frame
#'
#' @inherit plot_stacked_bar
#' @param sample_metadata a dataframe containing a SampleID column and additional columns for each property you want to add as metadata (data.frame)
#'
#' @return mutalisk dataframe with additional metadata columns (data.frame)
#' @export
#'
mutalisk_dataframe_add_metadata <- function(mutalisk_dataframe, sample_metadata){
  assertions::assert_names_include(mutalisk_dataframe, c("Signatures", "SampleID", "Contributions"))
  assertions::assert_names_include(sample_metadata, c("SampleID"))

  metadata_rows_orig=nrow(sample_metadata)
  sample_metadata <- sample_metadata %>%
    dplyr::distinct(.data[["SampleID"]], .keep_all = TRUE) %>%
    dplyr::filter(!is.na(.data[["SampleID"]]))

  if(metadata_rows_orig != nrow(sample_metadata))
    message("SampleID column in sample metadata table had duplicate SampleIDs or NA values. These have been removed")

  mutalisk_dataframe <- dplyr::left_join(x = mutalisk_dataframe, y = sample_metadata, by = "SampleID")
}

#' Mutalisk Dataframe
#'
#' Adds metadata from a file to the mutalisk dataframe
#'
#' @inherit plot_stacked_bar
#' @param metadata Either a path to csv file OR a dataframe. Must contain a header line which contains a SampleID column that matches that of mutalisk_dataframe (string)
#'
#' @return mutalisk dataframe with metadata columns (data.frame)
#' @export
mutalisk_dataframe_inform_user_of_metadata <- function(mutalisk_dataframe, metadata){

  if(is.character(metadata)){
    message("Metadata provided as a filepath ... checking if file exists")
    assertions::assert_file_exists(metadata)
    # Add assertion of read access

    message("    > Metadata file exists")
    metadata_df <- data.table::fread(file = metadata, header = TRUE)
  }
  else if(is.data.frame(metadata)){
    message("Metadata provided as a data.frame")
    metadata_df <- metadata
  }

  names(metadata_df) <- names(metadata_df) %>% sub(pattern = "^Tumor_Sample_Barcode$", replacement = "SampleID", x = .)
  checkmate::assert_names(colnames(metadata_df), must.include = "SampleID", .var.name = paste0("Metadata File Header: ", metadata))
  mutalisk_dataframe_add_metadata(mutalisk_dataframe, metadata_df)
}

#' Mutalisk Dataframe
#'
#' Get a vector of metadata columns from a mutalisk_dataframe.
#'
#' @inherit plot_stacked_bar
#'
#' @return a character vector containing names of metadata columns. If no metadata columns have been added, returns a zero length character vector. (character)
#' @export
#'
mutalisk_dataframe_metadata_column_names <- function(mutalisk_dataframe){
  checkmate::assert_names(names(mutalisk_dataframe), must.include = c("Signatures", "SampleID", "Contributions"))
  metadata_columns = colnames(mutalisk_dataframe) %>% purrr::keep(function(x) !(x %in% c("Signatures", "SampleID", "Contributions")))
  return(metadata_columns)
}

mutalisk_dataframe_metadata_column_message <- function(mutalisk_dataframe){
  col_names = mutalisk_dataframe_metadata_column_names(mutalisk_dataframe)
  message("Metadata columns: ", ifelse(length(col_names) ==0, "none. To add metadata, use the metadata argument of `mutalisk_best_signature_directory_to_dataframe`", paste(col_names, collapse = ", ")))
}



# Colors ------------------------------------------------------------------


colors=c("#84005e","#7ef631","#4913c1","#a7f915","#2d00a1","#01fa76","#ba58ff","#5ec400","#0053e8","#cae401","#9672ff","#beff53","#880095","#edff4d","#001f7a","#bcc700","#b4009a","#2affbd","#ff1f96","#00b150","#ff2b36","#00e0b1","#ff4c22","#29f5ff","#b90033","#8dffa2","#590062","#bdff8e","#320040","#ffd841","#0277e4","#fcae00","#54aeff","#7aa300","#d39cff","#018a20","#ff9ef8","#295500","#ff87c4","#003a00","#f0bfff","#002e0d","#ff9943","#06002d","#fdff9c","#2e0023","#b8ffb8","#a0004c","#00e8d3","#990035","#b1ffdb","#780012","#90e6ff","#932600","#01b1c1","#9e4f00","#87caff","#a27300","#015f9d","#ffe080","#003062","#e1ffbb","#001020","#fff9c7","#0b0f00","#ffeacf","#570009","#00a68a","#ff946d","#004348","#ff969c","#006c49","#abb4ff","#676900","#e1e1ff","#381900","#ffe5e5","#5f4900","#007d8d","#ffbab8")
short_color_set = c("8d2824","3c8dbd","ffd166","06d6a0","fffcf9")
#colors_v3_signatures = colors[1:77] %>% setNames()

#' Custom Palette
#'
#' @param n number of colors in pallete to return
#'
#' @return character vector of colors
#' @export
#'
#' @examples
#' pallette_cols23_customised(2)
pallette_cols23_customised <- function(n) {
  # derived from pals::cols25
  colors=c("#E84141","#1F78C8", "#778B00", "#6A33C2", "#ff7f00",
    "#FFD700", "#a6cee3", "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F", "#EEE685", "#C8308C", "#FF83FA", "#C814FA", "#0000FF",
    "#36648B", "#00E2E5", "#00FF00", "#33a02c", "#BEBE00", "#8B3B00",
    "#A52A3C")
  checkmate::assert_number(n, upper = length(colors), lower = 1)
  return(colors[1:n])
}

custom_pallete=function(n){
  return(colors[1:n])
}
