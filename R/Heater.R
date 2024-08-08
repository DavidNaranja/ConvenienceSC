#Wrapper for Seurat's AggregateExpression that outputs the fully processed matrix
processed_matrix <- function (obj, groupings = "ident", layer = "RNA"){
  mtx <- AggregateExpression(object = obj, group.by = groupings)[[layer]] %>% as.matrix

  zero_variance <- apply(mtx, 1, function(i) var(i) == 0)
  mtx <- mtx[!zero_variance, ]
  mtx <- na.omit(mtx)

  rsf <- mtx %>% t() %>% scale() %>% t()
  return (rsf)
}
