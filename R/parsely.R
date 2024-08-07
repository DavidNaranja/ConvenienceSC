parse_genes <- function(x){
  if(is.character(x)){
    if (length(x) == 0){
      character()
    }

    else{
      y <- (na.omit(unlist(stringr::str_split(x, ","))))
      y <- y[y != ""]
      return (y)
    }
  }
}

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
