ydx_draw_veen_plot <- function(data_list, picDir = paste0("veen_", GSE_number)) {
  num_vectors <- length(data_list)
  library(VennDiagram)
  if (num_vectors < 2 || num_vectors > 6) {
    stop("Input list must contain 2-6 vectors.")
  }
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  # Generate labels for the plot
  labels <- c(names(data_list))
  
  # Generate colors for the plot
  colors <- c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9")
  
  # Create a list of the input vectors
  input_vec_list <- lapply(data_list, as.character)
  
  # Create a Venn diagram using the appropriate method based on the number of input vectors
  if (num_vectors == 2) {
    venn.diagram(list(input_vec_list[[1]], input_vec_list[[2]]), fill = colors[1:2],
                 filename = file.path(picDir, "veen_plot.pdf"), main = "VEEN Plot",
                 cat.col = colors[1:2], category.names = labels[1:2], margin = 0.1)
  } else if (num_vectors == 3) {
    venn.diagram(list(input_vec_list[[1]], input_vec_list[[2]], input_vec_list[[3]]),
                 fill = colors[1:3], filename = file.path(picDir, "veen_plot.pdf"),
                 main = "VEEN Plot", cat.col = colors[1:3], category.names = labels[1:3],
                 margin = 0.1)
  } else if (num_vectors == 4) {
    venn.diagram(list(input_vec_list[[1]], input_vec_list[[2]], input_vec_list[[3]],
                      input_vec_list[[4]]), fill = colors[1:4],
                 filename = file.path(picDir, "veen_plot.pdf"), main = "VEEN Plot",
                 cat.col = colors[1:4], category.names = labels[1:4], margin = 0.1)
  } else if (num_vectors == 5) {
    venn.diagram(x=list(input_vec_list[[1]], input_vec_list[[2]], input_vec_list[[3]],
                        input_vec_list[[4]], input_vec_list[[5]]), fill = colors[1:5],
                 filename = file.path(picDir, "veen_plot.pdf"), main = "VEEN Plot",
                 cat.col = colors[1:5], category.names = labels[1:5], margin = 0.1)
  } else if (num_vectors == 6) {
    venn.diagram(list(input_vec_list[[1]], input_vec_list[[2]], input_vec_list[[3]],
                      input_vec_list[[4]], input_vec_list[[5]], input_vec_list[[6]]),
                 fill = colors[1:6], filename = file.path(picDir, "veen_plot.pdf"),
                 main = "VEEN Plot", cat.col = colors[1:6], category.names = labels[1:6],
                 margin = 0.1)
  }

  
  return()
}
