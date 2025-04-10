# ImageBinPlot

A Simple R function for binning and visualizing spatial transcriptomics data from Seurat objects. This utility supports both cell-based expression data and molecule-based transcript counts, with optional identity overlays using `scattermore`. Designed for Xenium data but adaptable to other spatial datasets.

## Features
- **Binning**: Aggregate spatial data into grids with customizable bin sizes.
- **Visualization**: Generate raster plots using `ggplot2` with the `viridis` color scale.
- **Data Types**: Supports two modes:
  - `expression`: For cell-level transcript counts.
  - `molecule`: For individual transcript coordinates.
- **Identity Overlay**: Optionally overlay cell identities (e.g., clusters) using `geom_scattermore`.
- **Multiple Fields of View**: Use `ImageBinPlotFOVs` to generate plots for multiple FOVs with optional same-scale adjustment.
- **Multi-Object Plotting**: Use `ImageBinPlotObjects` to visualize spatial data across multiple Seurat objects.
- **Customization**: Adjust bin size, color scaling, transparency, and layer order.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/simoncmo/ImageBinPlot.git
   ```

2. **Source the Script**:
   In your R session, source the script:
   ```R
   source("ImageBinPlot.R")
   ```

## Dependencies

Install the required R packages:
```R
install.packages(c("tidyverse", "Seurat", "ggplot2", "glue", "scattermore", "viridis", "patchwork"))
```
- `tidyverse`: For data manipulation and visualization.
- `Seurat`: For handling spatial transcriptomics data.
- `ggplot2`: For plotting.
- `glue`: For string interpolation.
- `scattermore`: For efficient scatter plot overlays.
- `viridis`: Provides color palettes.
- `patchwork`: For combining multiple plots.

## Usage

### Basic Setup
Load your Seurat object and source the script:
```R
obj <- readRDS("/path/to/your_seurat_object.rds")
library(qs)  # Optional if using .qs files
obj <- qread("/path/to/your_seurat_object.qs", nthreads = 30)
source("ImageBinPlot.R")
```

### Example 1: Expression-Based Binning
Plot binned expression data for a gene using `ImageBinPlot`:
```R
# Plot expression of "EPCAM" with 10µm bins
p <- ImageBinPlot(
  obj, 
  feature = "EPCAM", 
  bin_size = 10,  # in µm
  max_quantile = "q75", 
  type = "expression",
  # palette from viridis
  palette = 'viridis', # c('viridis', 'inferno', 'magma', 'plasma','cividis','mako','rocket','turbo')
  palette_begin = 0,
  palette_end = 1,
)
ggsave("EPCAM_expression.pdf", p, width = 10, height = 10)
```

### Example 2: Molecule-Based Binning
Plot binned transcript counts for a gene:
```R
# Plot molecule counts for "TP63" with 10µm bins
p <- ImageBinPlot(
  obj, 
  feature = "TP63", 
  bin_size = 10, 
  max_quantile = "q75", 
  type = "molecule"
)
ggsave("TP63_molecules.pdf", p, width = 10, height = 10)
```

### Example 3: With Identity Overlay
Overlay Seurat clusters on the expression plot:
```R
# Plot "EPCAM" expression with cluster overlay
p <- ImageBinPlot(
  obj, 
  feature = "EPCAM", 
  group.by = "seurat_clusters", 
  bin_size = 10, 
  max_quantile = "q75",
  type = "expression",
  ident_alpha = 0.3, 
  ident_pointsize = 1, 
  feature_on_top = TRUE
)
ggsave("EPCAM_with_clusters.pdf", p, width = 10, height = 10)
```

### Example 4: Filtered Identity Overlay
Show only specific clusters:
```R
# Overlay only cluster "0" on "EPCAM"
p <- ImageBinPlot(
  obj, 
  feature = "EPCAM", 
  group.by = "seurat_clusters", 
  filter_ident = "0", 
  bin_size = 10, 
  max_quantile = "q75", 
  type = "expression",
  ident_alpha = 0.8, 
  feature_on_top = TRUE
)
ggsave("EPCAM_cluster0.pdf", p, width = 10, height = 10)
```

### Example 5: Batch Processing for a Single Object
Generate plots for multiple genes:
```R
genes <- c("EPCAM", "TP63", "GeneX")
pdf("batch_plots.pdf", width = 10, height = 10)
for (gene in genes) {
  p <- ImageBinPlot(
    obj, 
    feature = gene, 
    bin_size = 10, 
    max_quantile = "q75", 
    type = "expression"
  )
  print(p)
}
dev.off()
```

### Example 6: Plotting Multiple Fields of View (FOVs)
Generate plots for several FOVs using `ImageBinPlotFOVs`:
```R
# Generate plots for all FOVs in the object with same scale
p_fovs <- ImageBinPlotFOVs(
  obj, 
  feature = "EPCAM", 
  bin_size = 10, 
  max_quantile = "q75",
  same_scale = TRUE
)
ggsave("EPCAM_all_FOVs.pdf", p_fovs, width = 10, height = 10)
```

### Example 7: Multi-Object Plotting
Combine plots for multiple spatial transcriptomics objects:
```R
# Assuming you have two Seurat objects: obj1 and obj2
p_objects <- ImageBinPlotObjects(
  obj_list = list(obj1 = obj1, obj2 = obj2), 
  feature = "EPCAM", 
  bin_size = 10, 
  max_quantile = "q75",
  same_scale = TRUE
)
ggsave("EPCAM_multi_object.pdf", p_objects, width = 10, height = 10)
```


## Notes
- Ensure your Seurat object(s) contain the required spatial data (centroids for expression or coordinate data for molecules).
- Adjust `max_quantile` to fine-tune color scaling.
- The script is designed with Xenium data in mind but can be adapted for other spatial data formats.

## Function Details

### Spatial Binning Plots
Visualize spatial transcriptomics data by creating binning plots for features (genes or molecules) in a Seurat object.

#### `ImageBinPlot`
Creates a spatial binning plot for a specified feature (gene or molecule) in a Seurat object.

```R
ImageBinPlot(
  obj,
  feature,
  group.by = NULL,
  x_col = "x",
  y_col = "y",
  fov = Images(obj)[[1]],
  x_min = NULL,
  x_max = NULL,
  y_min = NULL,
  y_max = NULL,
  bin_size = 10,
  max_multiplier = NULL,
  min = 0,
  max_quantile = "q75",
  type = "expression",
  assay = "Xenium",
  layer = "count",
  palette = c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket", "turbo"),
  palette_begin = 0,
  palette_end = 1,
  ident_alpha = 0.3,
  ident_pointsize = 1,
  filter_ident = NULL,
  feature_on_top = TRUE,
  flip_y = TRUE
)
```

##### Arguments
`obj`  
A Seurat object containing spatial transcriptomics data.

`feature`  
Character string specifying the gene or feature to visualize.

`group.by`  
Optional character string; column name in the Seurat object for identity overlay (e.g., cell type or cluster). Default is `NULL`.

`x_col`  
Character string specifying the column name for x-coordinate values. Default is `"x"`.

`y_col`  
Character string specifying the column name for y-coordinate values. Default is `"y"`.

`fov`  
Character string specifying the field of view to visualize; defaults to the first available FOV in the object (`Images(obj)[[1]]`).

`x_min`  
Numeric; minimum x-coordinate limit. If `NULL`, automatically set from the FOV’s centroid range.

`x_max`  
Numeric; maximum x-coordinate limit. If `NULL`, automatically set from the FOV’s centroid range.

`y_min`  
Numeric; minimum y-coordinate limit. If `NULL`, automatically set from the FOV’s centroid range.

`y_max`  
Numeric; maximum y-coordinate limit. If `NULL`, automatically set from the FOV’s centroid range.

`bin_size`  
Numeric; size of each spatial bin (in micrometers). Default is `10`.

`max_multiplier`  
Numeric; deprecated multiplier for color scale. Use `max_quantile` instead. Default is `NULL`.

`min`  
Numeric; minimum value for the color scale. Default is `0`.

`max_quantile`  
Character string (e.g., `"q75"`) specifying the quantile to determine the maximum value for the color scale from binned data. Default is `"q75"`.

`type`  
Character string specifying the type of data to plot; either `"expression"` (cell-based) or `"molecule"` (individual transcript counts). Default is `"expression"`.
When `type = "expression"` resulting average values per bin, while `type = "molecule"` values are count of molecules per bin.

`assay`  
Character string specifying the assay to use for data extraction from the Seurat object. Default is `"Xenium"`.

`layer`  
Character string specifying the data layer within the assay to use (e.g., `"count"`). Default is `"count"`.

`palette`  
Character vector of color palette options for visualization (e.g., `"viridis"`, `"inferno"`, etc.). Default is `c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket", "turbo")`.

`palette_begin`  
Numeric; start position of the selected color palette. Default is `0`.

`palette_end`  
Numeric; end position of the selected color palette. Default is `1`.

`ident_alpha`  
Numeric; transparency level for the identity overlay points, between 0 and 1. Default is `0.3`.

`ident_pointsize`  
Numeric; point size for the identity overlay. Default is `1`.

`filter_ident`  
Optional vector; filter to include only specific identities from `group.by`. Default is `NULL`.

`feature_on_top`  
Logical; if `TRUE`, the feature layer (binned data) is rendered on top of the identity overlay. Default is `TRUE`.

`flip_y`
Logical; Whether to flip the y-axis. Default is `TRUE`.

##### Value
A ggplot object representing the spatial binning plot.

---

#### `ImageBinPlotFOVs`
Generates bin plots for multiple fields of view (FOVs) in a Seurat object, with an option to enforce the same scale across FOVs.

```R
ImageBinPlotFOVs(
  obj,
  feature,
  group.by = NULL,
  fovs = NULL,
  bin_size = 10,
  min = 0,
  max_quantile = "q75",
  x_min = NULL,
  x_max = NULL,
  y_min = NULL,
  y_max = NULL,
  x_col = "x",
  y_col = "y",
  type = "expression",
  assay = "Xenium",
  layer = "count",
  ident_alpha = 0.3,
  ident_pointsize = 1,
  palette = c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket", "turbo"),
  palette_begin = 0,
  palette_end = 1,
  filter_ident = NULL,
  feature_on_top = TRUE,
  same_scale = FALSE,
  ncol = NULL
)
```

##### Arguments
`obj`  
A Seurat object containing spatial transcriptomics data.

`feature`  
Character string specifying the gene or feature to be visualized.

`group.by`  
Optional character string; column for identity overlay. Default is `NULL`.

`fovs`  
Vector of FOV names to plot. If `NULL`, defaults to all available FOVs in the object (`Images(obj)`).

`bin_size`  
Numeric; size of spatial bins (in micrometers). Default is `10`.

`min`  
Numeric; minimum value for the color scale. Default is `0`.

`max_quantile`  
Character string (e.g., `"q75"`) to determine the maximum value for the color scale. Default is `"q75"`.

`x_min`  
Numeric; optional minimum x-coordinate. Derived from data if `NULL`.

`x_max`  
Numeric; optional maximum x-coordinate. Derived from data if `NULL`.

`y_min`  
Numeric; optional minimum y-coordinate. Derived from data if `NULL`.

`y_max`  
Numeric; optional maximum y-coordinate. Derived from data if `NULL`.

`x_col`  
Character string specifying the column name for x-coordinate values. Default is `"x"`.

`y_col`  
Character string specifying the column name for y-coordinate values. Default is `"y"`.

`type`  
Character string; type of data to plot, either `"expression"` or `"molecule"`. Default is `"expression"`.

`assay`  
Character string; assay to use for data extraction. Default is `"Xenium"`.

`layer`  
Character string; data layer to use for extraction. Default is `"count"`.

`ident_alpha`  
Numeric; transparency level for identity overlay points, between 0 and 1. Default is `0.3`.

`ident_pointsize`  
Numeric; size of identity overlay points. Default is `1`.

`palette`  
Character vector of color palette options for plotting. Default is `c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket", "turbo")`.

`palette_begin`  
Numeric; start point for the color palette. Default is `0`.

`palette_end`  
Numeric; end point for the color palette. Default is `1`.

`filter_ident`  
Optional vector; filter to include only selected identities. Default is `NULL`.

`feature_on_top`  
Logical; if `TRUE`, the binned feature layer is drawn on top of the identity layer. Default is `TRUE`.

`same_scale`  
Logical; if `TRUE`, forces all FOV plots to use the same color scale. Default is `FALSE`.

`ncol`  
Numeric; optional number of columns when arranging multiple FOV plots in a grid. Default is `NULL`.

##### Value
A patchwork ggplot object combining bin plots for all specified FOVs.

---

#### `ImageBinPlotObjects`
Creates spatial binning plots for multiple Seurat objects, optionally across multiple FOVs, with an option for uniform scaling.

```R
ImageBinPlotObjects(
  obj_list,
  feature,
  group.by = NULL,
  fovs_list = NULL,
  bin_size = 10,
  min = 0,
  max_quantile = "q75",
  x_min = NULL,
  x_max = NULL,
  y_min = NULL,
  y_max = NULL,
  x_col = "x",
  y_col = "y",
  type = "expression",
  assay = "Xenium",
  layer = "count",
  ident_alpha = 0.3,
  ident_pointsize = 1,
  palette = c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket", "turbo"),
  palette_begin = 0,
  palette_end = 1,
  filter_ident = NULL,
  feature_on_top = TRUE,
  same_scale = FALSE,
  ncol = NULL
)
```

##### Arguments
`obj_list`  
A named list of Seurat objects to be plotted.

`feature`  
Character string; the gene or feature to visualize across objects.

`group.by`  
Optional character string; grouping variable for identity overlay across objects. Default is `NULL`.

`fovs_list`  
List of FOV vectors corresponding to each object. If `NULL`, defaults to all FOVs from each object via `Images(.x)`.

`bin_size`  
Numeric; size of the spatial bins (in micrometers). Default is `10`.

`min`  
Numeric; minimum value for the color scale. Default is `0`.

`max_quantile`  
Character string (e.g., `"q75"`) to determine the maximum value for the color scale. Default is `"q75"`.

`x_min`  
Numeric; optional minimum x-coordinate limit. Determined per object/FOV if `NULL`.

`x_max`  
Numeric; optional maximum x-coordinate limit. Determined per object/FOV if `NULL`.

`y_min`  
Numeric; optional minimum y-coordinate limit. Determined per object/FOV if `NULL`.

`y_max`  
Numeric; optional maximum y-coordinate limit. Determined per object/FOV if `NULL`.

`x_col`  
Character string; column name for x-coordinate values. Default is `"x"`.

`y_col`  
Character string; column name for y-coordinate values. Default is `"y"`.

`type`  
Character string; data type to visualize, either `"expression"` or `"molecule"`. Default is `"expression"`.

`assay`  
Character string; assay used for data extraction in each Seurat object. Default is `"Xenium"`.

`layer`  
Character string; data layer within the assay to use. Default is `"count"`.

`ident_alpha`  
Numeric; transparency level for identity overlay points, between 0 and 1. Default is `0.3`.

`ident_pointsize`  
Numeric; size of the identity overlay points. Default is `1`.

`palette`  
Character vector of color palette options for plotting. Default is `c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket", "turbo")`.

`palette_begin`  
Numeric; start point for the color palette. Default is `0`.

`palette_end`  
Numeric; end point for the color palette. Default is `1`.

`filter_ident`  
Optional vector; filter to restrict which identities to plot. Default is `NULL`.

`feature_on_top`  
Logical; if `TRUE`, the binned feature layer is drawn over the identity overlay. Default is `TRUE`.

`same_scale`  
Logical; if `TRUE`, applies a uniform color scale across all objects and their FOVs. Default is `FALSE`.

`ncol`  
Numeric; optional number of columns when arranging combined plots in a grid. Default is `NULL`.

##### Value
A patchwork ggplot object combining bin plots across all specified objects and FOVs.

--- 

This format aligns with the R documentation style you provided, including sections for function description, arguments, and return value, while maintaining clarity and consistency. Let me know if you'd like further adjustments!

## Contributing
Submit issues or pull requests to improve functionality or add features.

## License
MIT License - see LICENSE for details.
