
<!-- README.md is generated from README.Rmd. Please edit that file -->

# corplot

<!-- badges: start -->

<!-- badges: end -->

## Install

``` r
if(require("devtools")){
        devtools::install_github("cparsania/corplot")
} else{
        install.packages("devtools")
        devtools::install_github("cparsania/corplot")
}
```

## Correlation heatbox

## all vs all

``` r

expr_mat_file <- system.file("extdata" ,"example_data_expr_mat_01.txt" , package = "corplot")
expr_mat <- readr::read_delim(expr_mat_file , delim = "\t") 
cor_tbl <- corplot::get_pairwise_cor_tbl(expr_mat , var = "gene_name") 

cp <- corplot::get_corr_heat_box(cor_tbl,var1 = var1, var2 = var2 ,value = corr) 
cp + viridis::scale_fill_viridis() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Replicate A vs Replicate B

``` r
cor_tbl2 <- cor_tbl %>% dplyr::filter(grepl("Rep.A", var1) ) %>%  dplyr::filter(grepl("Rep.B", var2) )

corplot::get_corr_heat_box(cor_tbl2,var1 = var1, var2 = var2, value = corr) + 
  viridis::scale_fill_viridis() 
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## Pairwise scatter plot

``` r
groups_file <- expr_mat_file <- system.file("extdata" ,"example_data_01_sample_groups.txt" , package = "corplot")
groups <- readr::read_delim(file = groups_file,delim = "\t") 


csp <- corplot::get_pair_wise_scatter(dat_tbl = expr_mat, group_tbl = groups,var_plot = condition, var_plot_group = repl,dat_id = gene_name)

csp
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r


## Display corr value
cor_tbl2 <- cor_tbl %>% dplyr::rename(`Rep.A`=var1, `Rep.B` = var2) %>% 
  dplyr::filter(grepl("Rep.A" ,`Rep.A`)) %>%
  dplyr::filter(grepl("Rep.B" ,`Rep.B`)) %>% 
  TidyWrappers::tbl_replace_string("_.*" , "")

csp + ggplot2::geom_text(data = cor_tbl2,  x = 3, y = 18, ggplot2::aes(label = paste("r","=",corr , sep = "")) , 
                         fontface="italic" , col = "red")
```

<img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />
