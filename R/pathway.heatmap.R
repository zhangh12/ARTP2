

pathway.heatmap <- function(mat, row.names = rownames(mat), col.names = colnames(mat), 
                min.p = 1e-2, log.p = TRUE, na.color = 'grey90', reorder = TRUE, first.col.pathway = TRUE, 
                file = 'plot.png', width = 480, height = 480){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package reshape2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Package plyr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package grid needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if(log.p){
    mat <- -log10(mat)
  }
  
  id <- which(apply(mat, 2, function(x){!all(is.na(x))}))
  mat <- mat[, id, drop = FALSE]
  id <- which(apply(mat, 2, function(x){y <- x[!is.na(x)]; max(y) >= -log10(min.p)}))
  mat <- mat[, id, drop = FALSE]
  
  if(reorder){
    if(first.col.pathway){
      mat0 <- mat[, -1]
      tmp <- apply(mat0, 2, function(x){y <- x[!is.na(x)]; max(y)})
      mat0 <- mat0[, order(tmp, decreasing = TRUE), drop = FALSE]
      mat <- cbind(mat[, 1, drop = FALSE], mat0)
    }else{
      tmp <- apply(mat, 2, function(x){y <- x[!is.na(x)]; max(y)})
      mat <- mat[, order(tmp, decreasing = TRUE), drop = FALSE]
    }
  }
  
  rownames(mat) <- row.names
  colnames(mat) <- col.names
  
  nn <- ceiling(max(mat, na.rm = TRUE))
  
  mat <- as.data.frame(mat)
  mat$Name <- row.names
  mat$Name <- factor(mat$Name, levels = as.character(mat$Name))
  
  PTS <- 1:nrow(mat)
  
  mat$Name <- with(mat, reorder(Name, PTS))
  
  pos <- NULL
  for(i in 1:nrow(mat)){
    pos[i] <- which(levels(mat$Name) == mat$Name[i])
  }
  
  mat.m <- reshape2::melt(mat, id.vars = 'Name')
  mat.m <- plyr::ddply(mat.m, .(variable), transform)
  
  
  p1 <- ggplot(mat.m, aes(variable, Name)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradient(limits = c(1, nn), low = 'white', high = "red",
                        breaks = 1:nn, na.value = na.color)
  base_size <- 9
  p2 <- p1 + theme_grey(base_size = base_size) + 
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0), limits = levels(mat$Name)[pos]) + 
    theme(legend.position = "top", legend.key.size = unit(1, 'cm'), 
          legend.text = element_text(size = 10), legend.title = element_text(size = 12), 
          axis.ticks = element_blank(), 
          axis.text.x = element_text(size = base_size, 
                                     angle = 45, 
                                     hjust = 1, 
                                     vjust = 1, 
                                     colour = "black"), 
          axis.text.y = element_text(size = base_size, 
                                     colour = "black")) + 
    labs(fill = bquote(paste(-log[10], ' P-Value')))
  
  png(file, width = width, height = height)
  print(p2)
  dev.off()
  
}

