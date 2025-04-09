# spatial_binning_utils.R
# Utility functions for spatial binning and visualization of Seurat objects

# Required packages
library(tidyverse)
library(Seurat)
library(ggplot2)
library(glue)
library(scattermore)

# Helper function to reverse layer order in ggplot
reverse_layer_order <- function(plot_obj) {
  plot_obj$layers <- rev(plot_obj$layers)
  return(plot_obj)
}

# Get FOV centroid ranges
get_fov_centroid_ranges <- function(obj, fov = NULL) {
  fov <- fov %||% Images(obj)[[1]]
  fov_obj <- obj@images[[fov]]
  if (is.null(fov_obj)) {
    stop("FOV not found in the object.")
  }
  coord_df <- as.data.frame(fov_obj@boundaries$centroids@coords)
  list(
    x_min = min(coord_df$x),
    x_max = max(coord_df$x),
    y_min = min(coord_df$y),
    y_max = max(coord_df$y)
  )
}

# Fetch data with centroid coordinates
FetchDataWithCentroid <- function(obj, vars, assay = NULL, layer = "count", fov = NULL) {
  assay <- assay %||% DefaultAssay(obj)
  fov <- fov %||% Images(obj)[[1]]
  if (!assay %in% names(obj@assays)) {
    stop(paste("Assay", assay, "not found in the object."))
  }
  data <- FetchData(obj, vars = vars, assay = assay, layer = layer)
  coords <- obj@images[[fov]]@boundaries$centroids@coords
  as.data.frame(bind_cols(data, coords))
}

# Bin spatial data (supports both cell-based counts and molecule counts)
bin_spatial_data <- function(data, x_col = "x", y_col = "y", count_col = NULL, 
                             x_min = 0, x_max = 8000, y_min = 0, y_max = 8000, 
                             bin_size = 10, type = "expression") {
  
  # Create bin sequences
  x_bins <- seq(x_min, x_max, by = bin_size)
  y_bins <- seq(y_min, y_max, by = bin_size)
  
  # Assign points to bins
  data$x_bin <- cut(data[[x_col]], breaks = x_bins, include.lowest = TRUE, labels = FALSE)
  data$y_bin <- cut(data[[y_col]], breaks = y_bins, include.lowest = TRUE, labels = FALSE)
  
  if (type == "expression" && !is.null(count_col)) {
    # Sum counts for expression-based binning
    bin_df <- data %>%
      group_by(x_bin, y_bin) %>%
      summarise(count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
      filter(!is.na(x_bin) & !is.na(y_bin))
  } else if (type == "molecule") {
    # Count occurrences for molecule-based binning
    bin_counts <- table(data$x_bin, data$y_bin)
    bin_df <- as.data.frame.table(bin_counts)
    names(bin_df) <- c("x_bin", "y_bin", "count")
    
    # Ensure all bins are represented
    all_bins <- expand.grid(
      x_bin = seq_len(length(x_bins) - 1),
      y_bin = seq_len(length(y_bins) - 1)
    )
    bin_df <- merge(all_bins, bin_df, by = c("x_bin", "y_bin"), all.x = TRUE)
    bin_df$count[is.na(bin_df$count)] <- 0
  } else {
    stop("Invalid type specified. Use 'expression' or 'molecule'.")
  }
  
  # Convert bin indices to numeric and calculate coordinates
  bin_df$x_bin <- as.numeric(as.character(bin_df$x_bin))
  bin_df$y_bin <- as.numeric(as.character(bin_df$y_bin))
  bin_df$x_center <- x_min + (bin_df$x_bin - 0.5) * bin_size
  bin_df$y_center <- y_min + (bin_df$y_bin - 0.5) * bin_size
  
  # For molecule type, fill missing coordinates
  if (type == "molecule") {
    bin_df$x_center[is.na(bin_df$x_center)] <- x_min + (bin_df$x_bin[is.na(bin_df$x_center)] - 0.5) * bin_size
    bin_df$y_center[is.na(bin_df$y_center)] <- y_min + (bin_df$y_bin[is.na(bin_df$y_center)] - 0.5) * bin_size
  }
  
  return(bin_df)
}

# Main plotting function with optional identity overlay
# Old name: plot_binned_spatial_density 
# New name: ImageBinPlot
#' @param obj Seurat object
#' @param feature Gene to visualize
#' @param group.by Identity class to overlay
#' @param x_col Column name for x-coordinates
#' @param y_col Column name for y-coordinates
#' @param x_min Minimum x-coordinate
#' @param x_max Maximum x-coordinate
#' @param y_min Minimum y-coordinate
#' @param y_max Maximum y-coordinate
#' @param bin_size Size of bins in micrometers
#' @param max_multiplier Maximum multiplier for color scale
#' @param min Minimum value for color scale
#' @param type Type of data to visualize ("expression" or "molecule")
#' @param assay Assay to use for fetching data
#' @param layer Layer to use for fetching data
#' @param ident_alpha Alpha level for identity overlay
#' @param ident_pointsize Point size for identity overlay
#' @param filter_ident Optional filter for identity overlay
#' @param feature_on_top Whether to place feature layer on top of identity overlay
#' @return ggplot object
#' @export
#' @examples
#' ImageBinPlot(seurat_obj, feature = "GeneA", group.by = "cell_type", bin_size = 20)
#' ImageBinPlot(seurat_obj, feature = "GeneB", type = "molecule", bin_size = 50)
#' ImageBinPlot(seurat_obj, feature = "GeneC", group.by = "cell_type", filter_ident = c("Type1", "Type2"))
#' @description
#' This function creates a spatial binning plot for a specified gene or molecule in a Seurat object.
#' It allows for optional overlay of identity classes and customization of bin size, color scale, and other parameters.
#' #' The function uses ggplot2 for visualization and supports both expression and molecule data types.
#' #' @details
#' The function first checks if the specified gene exists in the object. It then fetches the data based on the specified type (expression or molecule) and bins the data into spatial bins of the specified size. The resulting plot shows the binned data with an optional overlay of identity classes, allowing for a clear visualization of spatial patterns in gene expression or molecule counts.
#' #' @return A ggplot object representing the spatial binning plot.
#' #' @export
#' ImageBinPlot
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_viridis_c scale_alpha_continuous theme_minimal labs coord_fixed theme
#' @importFrom dplyr mutate select filter group_by summarise ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom glue glue
#' @importFrom scattermore geom_scattermore
#' @importFrom scales squish
#' @importFrom purrr %||%
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_cols bind_rows arrange select filter group_by summarise 
#'   ungroup mutate distinct slice pull rename left_join right_join inner_join
#' @importFrom tidyr pivot_wider pivot_longer
ImageBinPlot <- function(obj, feature, group.by = NULL, 
                                        x_col = "x", y_col = "y", 
                                        x_min = NULL, x_max = NULL, 
                                        y_min = NULL, y_max = NULL, 
                                        bin_size = 10, max_multiplier = 0., min = 0, 
                                        type = "expression", assay = "Xenium", layer = "count",
                                        ident_alpha = 0.3, ident_pointsize = 1, 
                                        filter_ident = NULL, feature_on_top = TRUE) {
  
  # Check if the gene exists
  if (!feature %in% rownames(obj)) {
    warning(paste("Gene", feature, "not found in the object."))
    return(NULL)
  }
  
  # Fetch data based on type
  if (type == "expression") {
    input_data <- FetchDataWithCentroid(obj, vars = feature, assay = assay, layer = layer)
    count_col <- feature
  } else if (type == "molecule") {
    input_data <- obj@images$fov@molecules$molecules[[feature]] %>% as.data.frame
    count_col <- NULL
  } else {
    stop("Type must be 'expression' or 'molecule'.")
  }
  
  # Get coordinate ranges
  if (type == "expression") {
    range_df <- get_fov_centroid_ranges(obj)
  } else {
    range_df <- list(
      x_min = min(input_data$x),
      x_max = max(input_data$x),
      y_min = min(input_data$y),
      y_max = max(input_data$y)
    )
  }
  x_min <- x_min %||% range_df$x_min
  x_max <- x_max %||% range_df$x_max
  y_min <- y_min %||% range_df$y_min
  y_max <- y_max %||% range_df$y_max
  
  # Bin the data
  bin_data <- bin_spatial_data(
    input_data, x_col = x_col, y_col = y_col, count_col = count_col,
    x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max,
    bin_size = bin_size, type = type
  )
  
  # Base ggplot visualization
  p <- ggplot(bin_data, aes(x = x_center, y = y_center, fill = count, alpha = count)) +
    geom_raster() +
    scale_fill_viridis_c(limits = c(min, max(bin_data$count) * max_multiplier), oob = scales::squish) +
    scale_alpha_continuous(range = c(0.4, 1), limits = c(min, max(bin_data$count) * max_multiplier), oob = scales::squish) +
    theme_minimal() +
    labs(title = feature, subtitle = glue("({bin_size}µm x {bin_size}µm bins)")) +
    coord_fixed() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      text = element_text(color = "white"),
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 15)
    )
  
  # Add identity overlay if specified
  if (!is.null(group.by)) {
    df_ident <- FetchDataWithCentroid(obj, vars = group.by) %>%
      mutate(x_center = x, y_center = y) %>%
      select(-x, -y)
    
    if (!is.null(filter_ident)) {
      df_ident <- df_ident %>% filter(.data[[group.by]] %in% filter_ident)
    }
    
    p <- p + geom_scattermore(
      data = df_ident,
      inherit.aes = FALSE,
      aes(x = x_center, y = y_center, color = .data[[group.by]]),
      pointsize = ident_pointsize,
      alpha = ident_alpha
    )
    
    if (feature_on_top) {
      p <- reverse_layer_order(p)
    }
  }
  
  return(p)
}

# End of script