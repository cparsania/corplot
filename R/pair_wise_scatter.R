


#' Prepare pairwise scatter plot input data
#' @description Typical gene expression data stored in a tabular format, where rows are genes/features and columns are conditions/samples.
#' Except first column, which belongs to feature name, value in each cell correspond to expression magnitude of a respective feature.
#' Most gene expression data are produced in multiple replicate sets. Pairwise scatter plot is commonly used approach to check the quality of replicates.
#' Given the gene expression data as mentioned above along with sample groups, this function helps users to generate input data for pairwise scatter plot.
#' @param dat_tbl a tbl for which pairwise scatter plot to be generated. Rows congaing NA will be removed.
#' @param group_tbl a tbl containing sample groups. Refer details for more information on groups.
#' @param var_plot a variable name to be plotted.
#' @param var_plot_group a variable name to be used to group the variable \code{var_plot}.
#' @param dat_id a variable name from tbl \code{dat_tbl} storing feature ids. Typically name of first column form \code{dat_tbl}.
#' @importFrom rlang enquo quo sym quo_name
#' @importFrom glue glue
#' @importFrom dplyr pull left_join filter mutate select rename
#' @importFrom tidyr gather spread complete
#' @importFrom purrr map
#' @importFrom TidyWrappers tbl_remove_rows_NA_any
#' @importFrom ggplot2 ggplot geom_point aes facet_grid vars theme_bw theme element_text
#' @return a tbl, which can be used to generate pairwise scatter plot.
#' @details // TO DO
#'
get_pairwise_scatter_data <- function(dat_tbl , group_tbl , var_plot, var_plot_group, dat_id){

  # dat_tbl <- dat
  # group_tbl <- groups2
  # var_plot = rlang::quo(cond)
  # var_plot_group = rlang::quo(repl)
  # dat_id = rlang::quo(gene_name)

  var_plot <- rlang::enquo(var_plot)
  var_plot_group <- rlang::enquo(var_plot_group)
  dat_id <- rlang::enquo(dat_id)

  column_group_name  <-  rlang::quo(!!rlang::sym(colnames(group_tbl)[1]))

  ## input validations

  # 1) var_plot must present in group_tbl
  if(! rlang::quo_name(var_plot) %in% colnames(group_tbl)){
    stop(glue::glue("variable `{rlang::quo_name(var_plot)}` must present in tbl `group_tbl`") )
  }

  # 2) var_plot_group must present in group_tbl
  if(! rlang::quo_name(var_plot_group) %in% colnames(group_tbl)){
    stop(glue::glue("variable `{rlang::quo_name(var_plot_group)}` must present in tbl `group_tbl`") )
  }

  # 3) dat_id must present in dat_tbl
  if(! rlang::quo_name(dat_id) %in% colnames(dat_tbl)){
    stop(glue::glue("variable `{rlang::quo_name(dat_id)}` must present in tbl `dat_tbl`") )
  }

  # 4) candidates of group_tbl 1st column, must present in dat_tbl
  if(! group_tbl %>% dplyr::pull(1) %in% colnames(dat_tbl) %>% all()){
    stop(glue::glue("all candidates of `group_tbl` first column must present in `dat_tbl`") )
  }

  ## remove NA rows from input matrix
  dat_tbl <- dat_tbl %>% TidyWrappers::tbl_remove_rows_NA_any()
  group_tbl <- group_tbl %>% TidyWrappers::tbl_remove_rows_NA_any()

  ## long format
  dat_long <- dat_tbl  %>%
    TidyWrappers::tbl_remove_rows_NA_any() %>%
    tidyr::gather(vars, value , - !!dat_id)

  ## add groups
  dat_long_with_groups <- dat_long %>%
    dplyr::left_join(group_tbl , by = c("vars" = rlang::quo_name(column_group_name)))

  ## Remove NA from long data
  ## Logic: For a given column in dat_tbl if corresponding is not given group_tbl it will generate NA values in dat_long_with_groups matrix.
  ## Removing this NA make sure that only group candidates will be in final scatter plot

  dat_long_with_groups <- dat_long_with_groups  %>%
    TidyWrappers::tbl_remove_rows_NA_any()

  ## group candidates
  group_candidates <- dat_long_with_groups %>% dplyr::pull(!!var_plot_group) %>% unique()

  ## split by var_plot_group
  dat_by_var_group <- purrr::map(group_candidates, ~ dat_long_with_groups %>%
                                   dplyr::filter(!!var_plot_group == ..1) %>%
                                   dplyr::mutate(key = paste(!!dat_id , "_" ,!!var_plot,sep = "")) %>%
                                   dplyr::select(key,value)
  )
  names(dat_by_var_group) <- group_candidates

  plot_data <- dat_long_with_groups  %>%
    tidyr::spread(!!var_plot_group, value) %>%
    dplyr::mutate(combin_var = !!var_plot)  %>%
    tidyr::complete(!!dat_id, !!var_plot , combin_var)  %>%
    dplyr::mutate(key1 = paste(!!dat_id,"_",!!var_plot,sep="")) %>%
    dplyr::mutate(key2 = paste(!!dat_id,"_",combin_var,sep="")) %>%
    dplyr::select(!!dat_id,key1,key2, !!var_plot,combin_var) %>% unique()

  ## join and clean
  plot_data %>%
    dplyr::left_join(dat_by_var_group[[1]] , c("key1" = "key" ))  %>%
    dplyr::left_join(dat_by_var_group[[2]] , c("key2" = "key" )) %>%
    dplyr::select(-key1, -key2) %>%
    dplyr::rename( !!names(dat_by_var_group)[1] := !!var_plot) %>%
    dplyr::rename( !!names(dat_by_var_group)[2] := combin_var ) %>%
    dplyr::rename(!!paste(names(dat_by_var_group)[1],"value",sep="_") := "value.x") %>%
    dplyr::rename(!!paste(names(dat_by_var_group)[2] ,"value",sep="_") := "value.y")


}



#' Generate pairwise scatter plot
#' @description Generates pairwise scatter plot. One of the application of this function is to generates scatter plots between samples having 2 replicates.
#'
#' @param dat_tbl a tbl for which pairwise scatter plot to be generate.
#' @param group_tbl a tbl containing sample groups. Refer details for more information on groups.
#' @param var_plot a variable name, which is to be plotted.
#' @param var_plot_group a variable name, which is to be used to group the variable \code{var_plot}
#' @param dat_id a variable name from tbl \code{dat_tbl} storing feature ids. Typically name of first column form \code{dat_tbl}.
#' @param view_matrix logical, default TRUE, whether to display matrix view.
#'
#' @return an object of ggplot2
#' @export
#'
#'
get_pair_wise_scatter <- function(dat_tbl , group_tbl , var_plot, var_plot_group, dat_id , view_matrix =TRUE){

  group_tbl <- group_tbl
  dat_tbl <- dat_tbl
  var_plot <- rlang::enquo(var_plot)
  var_plot_group <- rlang::enquo(var_plot_group)
  dat_id <- rlang::enquo(dat_id)

  plot_data <-  get_pairwise_scatter_data(dat_tbl = dat_tbl , group_tbl = group_tbl , var_plot = !!var_plot , var_plot_group = !!var_plot_group, dat_id = !!dat_id)


  var_x <- rlang::sym(plot_data %>% colnames() %>% .[2])
  var_y <- rlang::sym(plot_data %>% colnames() %>% .[3])
  value_x <- rlang::sym(plot_data %>% colnames() %>% .[4])
  value_y <- rlang::sym(plot_data %>% colnames() %>% .[5])

  if(!view_matrix){
    plot_data <- plot_data %>% dplyr::filter(!! var_x == !!var_y)
  }

  # if(cor_val){
  #   cor_mat <- dat_tbl  %>%
  #     tidyr::gather(key, value, -!!dat_id) %>%
  #     dplyr::left_join(group_tbl , by = c("key" =  group_tbl %>% colnames()[1])) %>%
  #     dplyr::select(1,3,4,5) %>%
  #     tidyr::spread(!!var_plot,value) %>%
  #     dplyr::group_by(!!var_plot_group) %>%
  #     dplyr::summarise(corr = cor(Rep.A, Rep.B))
  # }

  gp <- plot_data %>%
    TidyWrappers::tbl_remove_rows_NA_any() %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = !!value_x , y = !!value_y )) +
    ggplot2:: theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20))

  if(view_matrix){
    gp <- gp + ggplot2::facet_grid(rows = ggplot2::vars(!!var_x),cols =  ggplot2::vars(!!var_y) )
  } else {
    gp <- gp + ggplot2::facet_wrap(c(ggplot2::vars(!!var_y)))
  }

}






