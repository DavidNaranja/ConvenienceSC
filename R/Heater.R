#Wrapper for Seurat's AggregateExpression that outputs the fully processed matrix
processed_matrix <- function (obj, groupings = "ident", layer = "RNA"){
  mtx <- Seurat::AggregateExpression(object = obj, group.by = groupings)[[layer]] %>% as.matrix

  zero_variance <- apply(mtx, 1, function(i) var(i) == 0)
  mtx <- mtx[!zero_variance, ]
  mtx <- na.omit(mtx)

  rsf <- mtx %>% t() %>% scale() %>% t()
  return(rsf)
}

#Wrapper for ComplexHeatmap's heatmap_annotation
#order of items in groupings vector will affect result
#For simplicity, aka until I confirm that all Seurat objects have name@meta.data
#ie. items are always treated as subclasses of the prior item
#What if length(groupings) == 0?
generate_heatmap_annotation <- function(obj_df, groupings, cols = NULL){
  column_df <- create_column_df(obj_df, groupings)

  cols_list <- create_cols_list(cols)

  ComplexHeatmap::HeatmapAnnotation(
    df = column_df,
    simple_anno_size = unit(3, "mm"),
    col = cols_list,
    gp = gpar(col = "white") ,
    annotation_name_gp = gpar(fontsize = 0),
    show_legend = T)
}

create_column_df <- function(obj_df, groupings){
  groups <- lapply(groupings, function(x){
    unique(obj_df[[x]])
  })
  num_unique_per_group <- unlist(lapply(groups, function(y){
    length(y)
  }))
  product_sum <- prod(num_unique_per_group)

  groups <- lapply(groups, function(x){
    rep(x, product_sum/length(x))
  })

  return (groups)
}
create_cols_list <- function(cols){
  if (is.null(cols)){
    lengths <- lapply(groupings, function(x){
      length(unique(obj_df[[x]]))
    })

    return (NULL) #Placeholder
  }
  else {
    return (cols)
  }

}

#Returns a function
fire_up_heatmap <- function(column_ann, rsf, col_ramp = colorRamp2(quantile(res_scaled, c(0.05, 0.5, 0.95)), c("blue", "black", "red"))){

  fire <- function(genes){
    setdiff(genes, rownames(rsf))
    res_scaled <- rsf[genes[genes %in% rownames(rsf)], ]

    gene_heatmap <- ComplexHeatmap::Heatmap(
      res_scaled, # Whether make cluster on columns or rows
      cluster_rows = T, #Cluster the genes according to the pre-defined gene dendrogram
      #row_km = 2,
      cluster_columns = F, #Cluster the sample according to the pre-defined sample dendrogram
      row_title = NULL, #Name the cluster
      row_title_rot = 0, #Rotate the row name
      top_annotation = column_ann,
      show_heatmap_legend = TRUE,
      show_row_names = T, row_names_side = "left",
      show_column_names = F,
      show_column_dend = F,
      column_title = t,
      col = col_ramp,
      column_dend_height = unit(0, "mm"),
      show_row_dend = F, #Not to put row dendrogram on the left
      row_names_gp = gpar(fontsize = 8, fontfamily = "sans", fontface = "bold"))
    draw(gene_heatmap)
  }
  return (fire)
}



