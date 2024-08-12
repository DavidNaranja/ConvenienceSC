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
generate_heatmap_annotation <- function(obj_df, groupings, cols = c("red", "blue", "purple", "yellow", "darkgreen", "red", "black", "lightgray")){
  gr <- lapply(groupings, function(x){
    unique(obj_df[[x]])
  })
  names(gr) <- groupings

  column_df <- create_column_df(gr)
  cols_list <- create_cols_list(gr, cols)

  ComplexHeatmap::HeatmapAnnotation(
    df = column_df,
    simple_anno_size = unit(3, "mm"),
    col = cols_list,
    gp = gpar(col = "white") ,
    annotation_name_gp = gpar(fontsize = 0),
    show_legend = T)
}

create_column_df <- function(groups){
  product_sum <- prod(lengths(groups)) #Gets the highest number of possible group configurations

  groups <- lapply(groups, function(x){
    rep(x, product_sum/length(x))
  })

  return (groups)
}
create_cols_list <- function(groups, cols){
  if (is.list(cols)){
    return (cols)
  }

  cols_col <- lapply(1:length(groups), function(i){
    group <- unlist(groups[i])
    group_length <- length(group)

    col_fun = colorRamp2::colorRamp2(c(0, 1), c(cols[i*2-1], cols[i*2]))
    group_cols <- lapply(1:group_length, function(y){
      col_fun(1/group_length * y)
    })
    group_cols <- unlist(group_cols)
    names(group_cols) <- group
    unlist(group_cols)
  })

  names(cols_col) <- names(groups)
  return (cols_col)
}

#Returns a function
fire_up_heatmap <- function(title, column_ann, rsf, colr =c("blue", "black", "red")){

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
      column_title = title,
      col = colorRamp2(quantile(res_scaled, c(0.05, 0.5, 0.95)), colr),
      column_dend_height = unit(0, "mm"),
      show_row_dend = F, #Not to put row dendrogram on the left
      row_names_gp = gpar(fontsize = 8, fontfamily = "sans", fontface = "bold"))
    draw(gene_heatmap)
  }
  return (fire)
}


