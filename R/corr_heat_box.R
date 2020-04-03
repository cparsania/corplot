



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
#' @importFrom ggplot2 ggplot geom_tile geom_label theme_bw theme
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
    ggplot2::geom_label(aes(x = !!var1, y = !!var2 , label = !!value)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 90))
  return(gp)

}

