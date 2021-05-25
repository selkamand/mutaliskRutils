
#' Mutalisk to dataframe
#'
#' You probably want to use tags$strong{mutaliskToDataFrame} instead.
#' See ?mutaliskToDataFrame
#'
#' @param mutalisk_files a vector of filepaths, each leading to the report.txt files output when downloading best_signature results for all vcfs in cohort
#'
#' @return tibble
#'
mutaliskToDataFrameSingleSample <- function(mutalisk_files){
  mutfile_v = readLines(mutalisk_files)
  sig_numbers_s = mutfile_v[grep(pattern = "\\*\\*THIS_SIGs ", x = mutfile_v)][1]
  sig_numbers_s = sig_numbers_s %>% strsplit(" ") %>% unlist()
  sig_numbers_f = sig_numbers_s[-1] %>% as.factor()

  sig_contributions_s = mutfile_v[grep(pattern = "\\*\\*THIS_SIG_CONTRIBUTIONS ", x = mutfile_v)][1]
  sig_contributions_s = sig_contributions_s %>% strsplit(" ") %>% unlist()
  sig_contributions_n = sig_contributions_s[-1] %>% as.numeric()

  patient_id_s = extractSampleNamesFromMutaliskFilenames(mutalisk_files)
  data.frame(SampleID = rep(patient_id_s, times = length(sig_numbers_f)), Signatures = sig_numbers_f, Contributions =  sig_contributions_n) %>%
    dplyr::tibble()
}

#' Mutalisk files to dataframe
#'
#' @param mutalisk_files a vector of filepaths, each leading to the report.txt files output when downloading best_signature results for all vcfs in cohort
#'
#' @return  a dataframe containing three columns:
#' \enumerate{
#' \item \strong{SampleID}: a sample identifier.
#' \item \strong{Signatures}: an identifier for a particular signature.
#' \item \strong{Contributions}: the percentage contribution of the signature to the patients genetic profile (0.1 = 10 percent).
#' }
#' @export
#'
mutaliskToDataFrame <- function(mutalisk_files){
  purrr::map_dfr(mutalisk_files,mutaliskToDataFrameSingleSample)
}

#' Sample Names From Mutalisk Output
#' @param mutalisk_filenames names of mutalisk output files (character)
#'
#' @return sample name (string)
#' @export
extractSampleNamesFromMutaliskFilenames  <- function(mutalisk_filenames){
  checkmate::assert_character(mutalisk_filenames)
  basename(mutalisk_filenames) %>%
    sub(pattern = ".*mutalisk_input_(.*?)\\..*$", replacement = "\\1", x =  mutalisk_filenames) %>%
    sub(pattern = "(.*)_.*$", replacement = "\\1", x = .) %>%
    return()
}

#' Mutalisk directory to dataframe
#'
#'
#' @param directory path to \strong{mutalisk_best_fit} folder.
#' To obtain, run your VCFs through mutalisk. Select Mutational Signature (Best only) and click 'Get the selected result for all samples at once'.
#' Then unzip the file, and youre ready to go
#'
#'
#' @return tibble
#' @export
#'
mutaliskBestSigatureDirectoryToDataframe <- function(directory){
 checkmate::assert_directory_exists(directory)

  filenames = dir(full.names = TRUE, path = directory, pattern = "mutalisk_input.*\\.txt$")
  mutaliskToDataFrame(filenames)
}


#' Cohort-Level Mutational Signature Visualisation
#'
#' A note of warning: for different mutalisk runs, this function will not enforce uniform colours for a single mutational signature. Better to get all data in at once and add a facet_wrap call
#'
#' @param mutalisk_dataframe a dataframe produced using mutaliskToDataFrame
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
#' @param facet_ncol number of columns in faceted plot (number)
#' @param fontsize_strip fontsize of facet titles (number)
#' @param fontsize_axis_title fontsize of axis titles (number)
#' @return a ggplot (gg)
#' @export
plotStackedBar <- function(mutalisk_dataframe, lump_type = "min_prop", lump_min=0.1, topn = 5, legend="right", legend_direction = NA, pal = pals::kovesi.diverging_rainbow_bgymr_45_85_c67, facet_column = NA, facet_ncol = 1, fontsize_strip = 18, fontsize_axis_title = 18){
  checkmate::assert_names(names(mutalisk_dataframe), must.include = c("Signatures", "SampleID", "Contributions"))
  checkmate::assert_choice(x = lump_type, choices = c("min_prop", "topn", "none"))
  checkmate::assert_choice(x = legend, choices = c("top", "left", "bottom", "right"))
  checkmate::assert_number(fontsize_strip)
  checkmate::assert_number(fontsize_axis_title)

  if(is.na(legend_direction))
    legend_direction = ifelse(legend=="right" | legend == "left", yes = "vertical", no = "horizontal")
  else
    checkmate::assert_choice(legend_direction, choices = c("horizontal", "vertical"))

  #Group Variables
  mutalisk_dataframe <- mutalisk_dataframe %>%
    dplyr::ungroup() %>%
    dplyr::group_by(SampleID)

  if (lump_type == "min_prop") {
    checkmate::assert_number(x = lump_min)
    message("Lumping together signatures with contributions < ", lump_min)
    mutalisk_dataframe <- mutalisk_dataframe %>%
      dplyr::group_by(SampleID) %>%
      dplyr::mutate(
        Signatures = forcats::fct_lump_min(Signatures, w = Contributions, min = 0.1)
      )
    mutalisk_dataframe$Signatures <- mutalisk_dataframe$Signatures %>%
      forcats::fct_relevel("Other", after = Inf) %>%
      droplevels()
  }

  if (lump_type == "topn") {
    checkmate::assert_number(x = topn, finite = TRUE)
    message(
      "Lumping together all signatures except those with the ",
      topn ,
      " greatest contributions"
    )
    mutalisk_dataframe <- mutalisk_dataframe %>%
      dplyr::group_by(SampleID) %>%
      dplyr::mutate(
        Signatures = forcats::fct_lump_n(Signatures, w = Contributions, n = topn)
      )
    mutalisk_dataframe$Signatures <- mutalisk_dataframe$Signatures %>%
      forcats::fct_relevel("Other", after = Inf) %>%
      droplevels()
  }

  if (lump_type == "none") {
    message("No signatures will be lumped together as 'other'. All signatures in a mutalisk 'best fit' set will be shown.")
  }

  gg = mutalisk_dataframe %>% ggplot2::ggplot() +
    ggplot2::geom_col(
      ggplot2::aes(
        x = SampleID,
        y = Contributions,
        fill = Signatures,
      ), position = "fill") +
        ggplot2::coord_flip() + ggplot2::xlab("Sample") + ggplot2::ylab("Contribution") +
        ggplot2::scale_fill_manual(
          values = c(
            pal(n = nlevels(mutalisk_dataframe$Signatures) - 1),"Grey"
          )
        ) +
        ggplot2::labs(fill = "Signature") +
        ggthemes::theme_fivethirtyeight() +
        ggplot2::theme(legend.position = legend, legend.direction = "vertical") +
        ggplot2::theme(axis.title = ggplot2::element_text(face = "bold", size = fontsize_axis_title)) +
        ggplot2::theme(strip.text= ggplot2::element_text(face = "bold", size = fontsize_strip), strip.background = ggplot2::element_blank()) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), plot.background = ggplot2::element_blank())

  if(!is.na(facet_column)){
    checkmate::assert(!checkmate::test_null(mutalisk_dataframe[[facet_column]]))
    checkmate::assert_number(facet_ncol)

    gg <- gg + ggplot2::facet_wrap(facets = facet_column, ncol = facet_ncol, scales = "free_y")
  }
      return(gg)
}

#' Cohort-Level Mutational Signature Visualisation
#'
#' @inherit plotStackedBar
#'
#' @return plotly visualisation
plotStackedBarInteractive <- function(mutalisk_dataframe, lump_type = "min_prop", lump_min=0.1, topn = 5, legend="top"){
  plotly::ggplotly(plotStackedBar(mutalisk_dataframe, lump_type = lump_type, lump_min=lump_min, topn = topn, legend=legend)) %>%
    return()
}


#' Plot Signature-Level Dotplot
#' Plots a signature-Level dotplot
#'
#' @inheritParams plotStackedBar
#'
#' @return a ggplot object
#' @export
#'
plotDistributionOfSignatureContributions <- function(mutalisk_dataframe){
  checkmate::assert_names(names(mutalisk_dataframe), must.include = c("Signatures", "SampleID", "Contributions"))

  levels(mutalisk_dataframe$Signatures) <- gtools::mixedsort(levels((mutalisk_dataframe$Signatures)))

  mutalisk_dataframe <- mutaliskDataframeExpand(mutalisk_dataframe)

  mutalisk_dataframe %>%
    ggplot2::ggplot(ggplot2::aes(Signatures, Contributions)) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::stat_summary(
      fun = median, fun.min = median, fun.max = median,
      geom = "crossbar", color = "red", width = 0.7, lwd = 0.4) +
    ggplot2::theme_bw()
}

#' Expand mutalisk_dataframe
#'
#' In normal mutalisk dataframe, each sample has data for ONLY the 1-7 signatures that comprise the 'best fit'.
#' This means that a signature level jitterplot won't show the samples where it was not included in this 'best fit' set, when we'd actually want to know that it contributed 0% to that sample.
#' This function fixes the issue by adding entries for ALL signature - sample pairs, with Contributions set to 0% where relevant.
#'
#' @inheritParams plotStackedBar
#'
#' @return dataframe containing all combinations of Sample ID and Signatures.
#' For cases where a signature was not included in the 'best fit' subset, Contribution is set to 0%.
#'
mutaliskDataframeExpand <- function(mutalisk_dataframe){
  checkmate::assert_names(names(mutalisk_dataframe), must.include = c("Signatures", "SampleID", "Contributions"))

  expand.grid(Signatures=unique(mutalisk_dataframe$Signatures),SampleID=unique(mutalisk_dataframe$SampleID)) %>%
      dplyr::mutate(Contributions = sapply(seq_along(Signatures),
                                           function(row_num) {
                                             res = mutalisk_dataframe[mutalisk_dataframe$Signatures == Signatures[row_num] & mutalisk_dataframe$SampleID == SampleID[row_num], "Contributions"] %>% unlist()
                                             ifelse(length(res) < 1, 0, res)
                                             }))
}

