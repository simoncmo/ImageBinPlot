# spatial_binning_utils.R
# Utility functions for spatial binning and visualization of Seurat objects

# Required packages
library(tidyverse)
library(Seurat)
library(ggplot2)
library(glue)
library(scattermore)
library(viridis)
library(patchwork)

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


# Helper function to get cut values from mulitple plots
GetNumericQuantile <- function(quantile_string){
	str_remove(quantile_string, 'q') %>% as.numeric() / 100
}

GetPlotValueQuantile <- function(p_list, value_column = 'count', quantile = 0.75, remove_zero=TRUE) {
	# Extract the data from each plot
	values <- map(p_list, ~.$data[[value_column]]) %>% reduce(c)
	# Remove 0 
	if(remove_zero) values <- values[values > 0]
	
	return(as.numeric(quantile(values, probs = quantile)))
}

# Helper function for getting viridis palette
viridisColor <- function(
  palette = c('viridis', 'inferno', 'magma', 'plasma','cividis','mako','rocket','turbo'),
  begin = 0.3,
  end = 1,
  n = 10,
  alpha = 1,
  ...
) {
  palette <- match.arg(palette)
  get(palette)(n = n, begin = begin, end = end, alpha = alpha, ...)
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
#' @param max_multiplier Maximum multiplier for color scale - Deprecated
#' @param max_quantile Maximum quantile for color scale
#' @param fov Field of view to visualize
#' @param palette Color palette to use
#' @param palette_begin Start of color palette
#' @param palette_end End of color palette
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
                                        fov = NULL,
                                        x_min = NULL, x_max = NULL, 
                                        y_min = NULL, y_max = NULL, 
                                        bin_size = 10, 
                                        max_multiplier = NULL, 
                                        min = 0, 
                                        max_quantile = "q75",
                                        type = "expression", assay = "Xenium", layer = "count",
                                        # palette
                                        palette = c('viridis', 'inferno', 'magma', 'plasma','cividis','mako','rocket','turbo'),
                                        palette_begin = 0,
                                        palette_end = 1,
                                        ident_alpha = 0.3, ident_pointsize = 1, 
                                        filter_ident = NULL, feature_on_top = TRUE) {
  
  # Check if the gene exists
  if (!feature %in% rownames(obj)) {
    warning(paste("Gene", feature, "not found in the object."))
    return(NULL)
  }

  # get Fov
  fov <- fov %||% Images(obj)[[1]]
  if (!fov %in% Images(obj)) {
    warnings(paste("FOV", fov, "not found in the object."))
    return(NULL)
  }
  
  # Fetch data based on type
  if (type == "expression") {
    input_data <- FetchDataWithCentroid(obj, vars = feature, assay = assay, layer = layer)
    count_col <- feature
  } else if (type == "molecule") {
    input_data <- obj@images[[fov]]@molecules$molecules[[feature]] %>% as.data.frame
    count_col <- NULL
  } else {
    stop("Type must be 'expression' or 'molecule'.")
  }
  
  # Get coordinate ranges
  if (type == "expression") {
    range_df <- get_fov_centroid_ranges(obj, fov = fov)
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
  
# Get the quantile value for the max multiplier
  max_quantile_numeric <- GetNumericQuantile(max_quantile)
  max_value <- quantile(bin_data$count %>% .[.>0], probs = max_quantile_numeric)
  print(max_quantile_numeric)
  print(max_value)
  # Set the maximum value for the color scale
  if (!is.null(max_multiplier)) {
    warning("Deprecation warning: max_multiplier is deprecated. Use max_quantile = 'q75' instead.")
  }
  
  # Base ggplot visualization
  p <- ggplot(bin_data, aes(x = x_center, y = y_center, fill = count, alpha = count)) +
    geom_raster() +
    scale_fill_gradientn(limits = c(min, max_value), colors = viridisColor(palette, n = 10, begin = palette_begin, end = palette_end), oob = scales::squish) +
    scale_alpha_continuous(range = c(0.4, 1), limits = c(min, max_quantile_numeric), oob = scales::squish) +
    theme_minimal() +
    labs(title = feature, subtitle = fov, caption = glue("({bin_size}µm x {bin_size}µm bins)")) +
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
    ) +
    # set legend title to either BinnedExpression or BinnedMolecule
    labs(fill = ifelse(type == "expression", "BinnedExpression", "BinnedMolecule")) +
    labs(alpha = ifelse(type == "expression", "BinnedExpression", "BinnedMolecule")) 
  
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


# helper function to seperate combined plot to a list of plot
CombinedPlotToList <- function(combined_plot) {
  # Extract the individual plots from the combined plot
  n_plots <- length(combined_plot)
  plot_list <- vector("list", n_plots)
  names(plot_list) <- paste0("Plot", seq_len(n_plots))
  for (i in seq_len(n_plots)) {
    plot_list[[i]] <- combined_plot[[i]]
  }
  
  # Return the list of plots
  return(plot_list)
}

# Helper functino to Make Plots same scale
MakeCombinePlotSameScale <- function(combined_plot, min = 0, max_quantile = "q75", colors = NULL) {
  # Get the quantile value for the max multiplier
  max_quantile_numeric <- GetNumericQuantile(max_quantile)
  max_value <- GetPlotValueQuantile(CombinedPlotToList(combined_plot), value_column = 'count', quantile = max_quantile_numeric, remove_zero=TRUE)
  
  # Set color
  colors <- colors %||% viridisColor('viridis', n = 10, begin = 0, end = 1) # default viridis color
  
  # Set the same scale for all plots
  combined_plot <- combined_plot & 
    scale_fill_gradientn(limits = c(min, max_value), colors = colors, oob = scales::squish) &
    scale_alpha_continuous(range = c(0.4, 1), limits = c(min, max_quantile_numeric), oob = scales::squish) 
  
  return(combined_plot)
}


#' ImageBinPlotFOVs
#'
#' This function generates bin plots for multiple fields of view (FOVs) from a spatial transcriptomics object. 
#' If no FOVs are specified, it defaults to using all FOVs available in the object. The function allows customization 
#' of bin size, color palette, and scaling options.
#'
#' @param obj A spatial transcriptomics object containing the data to be plotted.
#' @param feature The feature to be visualized (e.g., gene expression).
#' @param group.by Optional grouping variable for coloring points.
#' @param fovs A vector of FOVs to be plotted. Defaults to NULL, which uses all FOVs in the object.
#' @param bin_size The size of the bins for aggregating data. Default is 10.
#' @param min Minimum value for the color scale. Default is 0.
#' @param max_quantile Maximum quantile for the color scale. Default is "q75".
#' @param x_min, x_max, y_min, y_max Optional limits for the x and y axes.
#' @param x_col, y_col Column names for x and y coordinates. Default are "x" and "y".
#' @param type The type of data to plot (e.g., "expression"). Default is "expression".
#' @param assay The assay to use for data extraction. Default is "Xenium".
#' @param layer The data layer to use. Default is "count".
#' @param ident_alpha Transparency level for points. Default is 0.3.
#' @param ident_pointsize Size of points in the plot. Default is 1.
#' @param palette A vector of color palettes to choose from. Default includes 'viridis', 'inferno', etc.
#' @param palette_begin, palette_end Range of the color palette to use. Default is 0 to 1.
#' @param filter_ident Optional filter for specific identities.
#' @param feature_on_top Logical, whether to plot the feature layer on top. Default is TRUE.
#' @param same_scale Logical, whether to use the same scale for all FOVs. Default is FALSE.
#' @param ncol Number of columns for arranging plots. Default is NULL.
#'
#' @return A combined plot of all FOVs, optionally with the same scale.
#' @export
#'
#' @examples
#' # Generate plots for all FOVs in an object
#' ImageBinPlotFOVs(obj = my_object, feature = "GeneA")
#'
#' # Generate plots for specific FOVs with custom bin size
#' ImageBinPlotFOVs(obj = my_object, feature = "GeneA", fovs = c("FOV1", "FOV2"), bin_size = 20)
#'
#' # Generate plots with the same scale across all FOVs
#' ImageBinPlotFOVs(obj = my_object, feature = "GeneA", same_scale = TRUE)
#'
#' @seealso \code{\link{ImageBinPlot}}
#'
#' @importFrom patchwork wrap_plots
#' @importFrom purrr map
#'



# Same scale for All fov provided
# call ImageBinPlot with fovs = NULL. default to all fovs using Images(obj) to get all fov
ImageBinPlotFOVs <- function(obj, feature, group.by = NULL, 
                               fovs = NULL,
                               bin_size = 10, min = 0, max_quantile = "q75", 
                               x_min = NULL, x_max = NULL, 
                               y_min = NULL, y_max = NULL, 
                               x_col = "x", y_col = "y", 
                               type = "expression", assay = "Xenium", layer = "count",
                               ident_alpha = 0.3, ident_pointsize = 1, 
                               # palette
                               palette = c('viridis', 'inferno', 'magma', 'plasma','cividis','mako','rocket','turbo'),
                               palette_begin = 0,
                               palette_end = 1,
                               filter_ident = NULL, feature_on_top = TRUE,
                               same_scale = FALSE, ncol = NULL) {
  
  # Get all FOVs if none specified
  if (is.null(fovs)) {
    fovs <- Images(obj)
  }
  
  # Create plots for each FOV
  plot_list <- map(fovs, function(fov){
    p <- ImageBinPlot(
      obj = obj, feature = feature, group.by = group.by,
      x_col = x_col, y_col = y_col, fov = fov,
      x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max,
      bin_size = bin_size, 
      min = min,
      max_quantile = max_quantile,
      type = type, assay = assay, layer = layer,
      ident_alpha = ident_alpha, ident_pointsize = ident_pointsize,
      filter_ident = filter_ident, feature_on_top = feature_on_top
    )

    return(p)
  })

  # Combine plots
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol)
  
  # Make Same scale 
  if (same_scale) {combined_plot <- MakeCombinePlotSameScale(combined_plot, min = min, max_quantile = max_quantile, 
    colors = viridisColor(palette, n = 10, begin = palette_begin, end = palette_end) 
    ) }
  

  return(combined_plot)
}


#' ImageBinPlotObjects
#'
#' This function generates bin plots for multiple spatial transcriptomics objects, optionally across multiple FOVs. 
#' It allows customization of bin size, color palette, and scaling options. If no FOVs are specified, it defaults to 
#' using all FOVs available in each object.
#'
#' @param obj_list A named list of spatial transcriptomics objects to be plotted.
#' @param feature The feature to be visualized (e.g., gene expression).
#' @param group.by Optional grouping variable for coloring points.
#' @param fovs_list A list of FOVs for each object. Defaults to NULL, which uses all FOVs in each object.
#' @param bin_size The size of the bins for aggregating data. Default is 10.
#' @param min Minimum value for the color scale. Default is 0.
#' @param max_quantile Maximum quantile for the color scale. Default is "q75".
#' @param x_min, x_max, y_min, y_max Optional limits for the x and y axes.
#' @param x_col, y_col Column names for x and y coordinates. Default are "x" and "y".
#' @param type The type of data to plot (e.g., "expression"). Default is "expression".
#' @param assay The assay to use for data extraction. Default is "Xenium".
#' @param layer The data layer to use. Default is "count".
#' @param ident_alpha Transparency level for points. Default is 0.3.
#' @param ident_pointsize Size of points in the plot. Default is 1.
#' @param palette A vector of color palettes to choose from. Default includes 'viridis', 'inferno', etc.
#' @param palette_begin, palette_end Range of the color palette to use. Default is 0 to 1.
#' @param filter_ident Optional filter for specific identities.
#' @param feature_on_top Logical, whether to plot the feature layer on top. Default is TRUE.
#' @param same_scale Logical, whether to use the same scale for all objects. Default is FALSE.
#' @param ncol Number of columns for arranging plots. Default is NULL.
#'
#' @return A combined plot of all objects and their respective FOVs, optionally with the same scale.
#' @export
#'
#' @examples
#' # Generate plots for all objects and their FOVs
#' ImageBinPlotObjects(obj_list = list(obj1 = my_object1, obj2 = my_object2), feature = "GeneA")
#'
#' # Generate plots for specific FOVs in each object
#' ImageBinPlotObjects(obj_list = list(obj1 = my_object1, obj2 = my_object2), 
#'                     feature = "GeneA", 
#'                     fovs_list = list(c("FOV1", "FOV2"), c("FOV3", "FOV4")))
#'
#' # Generate plots with the same scale across all objects
#' ImageBinPlotObjects(obj_list = list(obj1 = my_object1, obj2 = my_object2), 
#'                     feature = "GeneA", 
#'                     same_scale = TRUE)
#'
#' @seealso \code{\link{ImageBinPlotFOVs}}
#'
#' @importFrom patchwork wrap_plots
#' @importFrom purrr pmap
#' @importFrom purrr reduce
#'
ImageBinPlotObjects <- function(
  obj_list, feature, group.by = NULL, 
  fovs_list = NULL,
  bin_size = 10, min = 0, max_quantile = "q75", 
  x_min = NULL, x_max = NULL, 
  y_min = NULL, y_max = NULL, 
  x_col = "x", y_col = "y", 
  type = "expression", assay = "Xenium", layer = "count",
  ident_alpha = 0.3, ident_pointsize = 1, 
  # palette
  palette = c('viridis', 'inferno', 'magma', 'plasma','cividis','mako','rocket','turbo'),
  palette_begin = 0,
  palette_end = 1,
  filter_ident = NULL, feature_on_top = TRUE,
  same_scale = FALSE, ncol = NULL) {
  # Set object name if not provided
  if (is.null(names(obj_list))) {
    names(obj_list) <- paste0("Object", seq_along(obj_list))
  }
  p_list_all <- pmap(
    list(
      object_name = names(obj_list),
      object_use = obj_list,
      fovs_use = fovs_list %||% map(obj_list, ~ Images(.x)) # get all fovs from each object
    ),function(object_name, object_use, fovs_use){
        p_combined <- ImageBinPlotFOVs(
          obj = object_use, feature = feature, group.by = group.by,
          fovs = fovs_use,
          bin_size = bin_size, min = min, max_quantile = max_quantile, 
          x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max, 
          x_col = x_col, y_col = y_col, 
          type = type, assay = assay, layer = layer,
          ident_alpha = ident_alpha, ident_pointsize = ident_pointsize,
          filter_ident = filter_ident, feature_on_top = feature_on_top,
          same_scale = FALSE, # Same scale will be done at all object level
          ncol = ncol
        )
        # Turn combined plot to plot_list
        return(CombinedPlotToList(p_combined))
    }
  ) %>% reduce(c) # combine all plot_list to a single list
  
  # Combine plots
  combined_plot <- patchwork::wrap_plots(p_list_all, ncol = ncol)
  
  # Make Same scale 
  if (same_scale) {combined_plot <- MakeCombinePlotSameScale(combined_plot, min = min, max_quantile = max_quantile, 
          colors = viridisColor(palette, n = 10, begin = palette_begin, end = palette_end)
    ) 
    }

  return(combined_plot)
  
  }