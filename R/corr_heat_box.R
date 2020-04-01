## Corr heatbox between sample pairs. samples must not repeat on any axis



# dat <- tibble::tribble(
#   ~genes, ~cg_xbp1_thp1_0h_polII_set1, ~cg_xbp1_thp1_2h_polII_set1, ~cg_xbp1_thp1_4h_polII_set1, ~cg_xbp1_thp1_6h_polII_set1, ~cg_xbp1_thp1_8h_polII_set1, ~cg_xbp1_thp1_0h_polII_set2, ~cg_xbp1_thp1_2h_polII_set2, ~cg_xbp1_thp1_4h_polII_set2, ~cg_xbp1_thp1_6h_polII_set2, ~cg_xbp1_thp1_8h_polII_set2,
#   "CAGL0A00105g",                 5.590530589,                 6.243382014,                 6.749620899,                 7.825111755,                 6.100145757,                 4.174793084,                  10.1650902,                 5.481785315,                 4.770598636,                 6.965798885,
#   "CAGL0A00110g",                  12.8425952,                 18.20034983,                 21.81659112,                 24.64828544,                 22.75378761,                 12.49357126,                 21.87086403,                 22.11528107,                 16.57443902,                 24.65169314,
#   "CAGL0A00116g",                 16.30495365,                 23.75513046,                 28.96951606,                 32.50976156,                 30.71185412,                 16.44640396,                 27.37725284,                 29.98144642,                 22.22565767,                 33.10066498,
#   "CAGL0A00132g",                 4.637555142,                 5.409842409,                 5.961671818,                 7.792807562,                 5.484018909,                 3.124362835,                 6.586119741,                 8.801859174,                 4.785924245,                 3.313132577,
#   "CAGL0A00154g",                 7.964959452,                 4.557837917,                 6.438378282,                 4.882072028,                 5.389297234,                 10.55931086,                 8.605121592,                 4.559572162,                    4.297054,                 1.986698259,
#   "CAGL0A00165g",                 8.868599914,                 7.537428691,                 9.537553437,                 8.350892851,                 9.118221384,                 10.19861793,                 7.937707141,                 10.97718701,                 8.644499961,                 5.756728672
# )
#



#' Prepare pairwise corr tibble
#' @description In general gene expression data are stored in tabular format, where first column is gene / feature name while rest other are expression values in respective samples/conditions.
#' Given such a data this function returns a long format tbl containing correlation value for each pair of variable.
#' @param dat a tbl.
#' @param var a column name denoting rownames. Default "genes".
#' @param transform logical, default TRUE, whether to log transform (log2) numeric columns.
#' @param method same as argument \code{method} in \code{stats::cor}.
#'
#' @return a tbl.
#' @export
#' @importFrom dplyr mutate_if
#' @importFrom TidyWrappers tbl_remove_rows_NA_any
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble
#' @importFrom tidyr gather
#' @importFrom rlang enquo
#' @importFrom ggplot2 ggplot geom_tile geom_label
#'
get_pairwise_cor_tbl <- function(dat , var = "genes", transform = TRUE, method = "pearson"){


  ## transform
  if(transform){
    dat <- dat %>% dplyr::mutate_if(is.numeric , ~(log2(. + 1)))
  }

  ## pair wise corr matrix
  cor_mat <- dat %>%
    TidyWrappers::tbl_remove_rows_NA_any() %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var) %>%
    as.matrix() %>%
    cor(method = method) %>%
    round(2)

  ## clean tibble
  cor_tbl <- cor_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var ="var1") %>%
    tibble::as_tibble() %>%
    tidyr::gather(key = "var2" , value= "corr" ,  -var1)

  return(cor_tbl)
}



#' Generate pairwise correlation heatbox
#' @description Given a tibble of pairwise correlation get a correlation heatbox.
#' @param pairwise_cor_tbl a tbl containing three columns
#' \enumerate{
#' \item \code{column 1)} values to show on x-axis.
#' \item \code{column 2)} values to show on y-axis.
#' \item \code{column 2)} values denoting correlation between \code{column 1} and \code{column 2}.
#' }
#' @param var1 a column to show on x-axis.
#' @param var2 a column to show on y-axis.
#' @param value a column denoting corr values.
#'
#' @return an object of ggplot
#' @export
#'
get_corr_heat_box <- function(pairwise_cor_tbl , var1, var2, value ){

  var1  = rlang::enquo(var1)
  var2  = rlang::enquo(var2)
  value = rlang::enquo(value)

  ## corr heat box
  gp <- pairwise_cor_tbl %>%
    ggplot2::ggplot() +
    ggplot2::geom_tile(aes(x = !!var1, y = !!var2 , fill = !!value)) +
    ggplot2::geom_label(aes(x = !!var1, y = !!var2 , label = !!value))
  return(gp)

}

