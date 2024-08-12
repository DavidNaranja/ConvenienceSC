#' Get list of genes
#'
#' @param x character containing genes to grab
#' @param keep character list of genes to keep
#' @param sep seperator of gene names
#'
#' @return List of genes. Only those in keep if keep is not NULL
#' @export
#'
#' @description
#' Takes a single string of gene names and turns them into a list
#'
#'
#' @examples
#' parse_genes("Xist,cdk9,Est3, Pea3")
parse_genes <- function(x, keep = NULL, sep = ","){
  if(is.character(x)){
    if (length(x) == 0){
      character()
    }
    y <- (na.omit(unlist(stringr::str_split(x, sep))))
    y <- y[y != ""]

    if (is.character(keep)){
      y <- y[y %in% keep]
    }

    return (y)
  }
}



#' Gets list of genes from gsea tsv
#'
#' @param x String path to tsv
#'
#' @return List of genes from tsv
#' @export
#'
#' @examples
#' #parse_gsea_tsv("../data/BIOCARTA_ETS_PATHWAY.v2023.2.Mm.tsv") Get more complete examples
parse_gsea_tsv <- function(x){
  stopifnot(is.character(x))

  g <- readr::read_tsv(x)
  r <- colnames(g)[2]
  g <- g[[r]][17]
  return (parse_genes(g))
}

parse_rds_vector <- function(x){
  parse_genes(readRDS(x))
}


#' Grab all genes from a url text file
#'
#' @param x String url
#' @param re Regex expression to use, default looks for the text ".symbol." preceding the gene symbol
#'
#' @return A character vector of all found genes
#' @export
#'
#' @examples
#' #Fill later
parse_url_raw <- function(x, re = "(?<=(\\.symbol\\.))[:alnum:]+"){
  y <- url(x) %>%
    readLines() %>%
    textConnection() %>%
    read.table(header = TRUE, sep = ',') %>%
    colnames() %>%
    str_extract(re) %>%
    na.omit() %>%
    str_to_title()
  return (y[y != ""])
}
