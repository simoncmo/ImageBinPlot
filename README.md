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

### `ImageBinPlot`
Creates a spatial binning plot for a specified feature (gene or molecule) in a Seurat object.  
| Parameter         | Type       | Default Value                              | Description                                                                                                                |
|-------------------|------------|--------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| obj               | Seurat     | —                                          | A Seurat object containing spatial transcriptomics data.                                                                |
| feature           | character  | —                                          | The gene or feature to visualize.                                                                                        |
| group.by          | character  | NULL                                       | Optional; column name in the Seurat object for identity overlay (e.g. cell type or cluster).                               |
| x_col             | character  | "x"                                        | Column name for x-coordinate values.                                                                                      |
| y_col             | character  | "y"                                        | Column name for y-coordinate values.                                                                                      |
| fov               | character  | Images(obj)[[1]]                           | Field of view to visualize; defaults to the first available FOV in the object.                                           |
| x_min             | numeric    | NULL                                       | Minimum x-coordinate limit; if NULL, automatically set from the FOV’s centroid range.                                    |
| x_max             | numeric    | NULL                                       | Maximum x-coordinate limit; if NULL, automatically set from the FOV’s centroid range.                                    |
| y_min             | numeric    | NULL                                       | Minimum y-coordinate limit; if NULL, automatically set from the FOV’s centroid range.                                    |
| y_max             | numeric    | NULL                                       | Maximum y-coordinate limit; if NULL, automatically set from the FOV’s centroid range.                                    |
| bin_size          | numeric    | 10                                         | Size of each spatial bin (in micrometers).                                                                               |
| max_multiplier    | numeric    | NULL                                       | (Deprecated) Previously used multiplier for color scale; use max_quantile instead.                                       |
| min               | numeric    | 0                                          | Minimum value for the color scale.                                                                                       |
| max_quantile      | character  | "q75"                                      | Quantile string (e.g., "q75") to determine the maximum value for the color scale from the binned data.                     |
| type              | character  | "expression"                               | Type of data to plot; either "expression" (cell-based) or "molecule" (individual transcript counts).                     |
| assay             | character  | "Xenium"                                   | Assay to use for data extraction from the Seurat object.                                                                 |
| layer             | character  | "count"                                    | Data layer within the assay to use (commonly "count").                                                                   |
| palette           | character  | c('viridis','inferno','magma','plasma','cividis','mako','rocket','turbo') | Color palette to use for visualization.                                 |
| palette_begin     | numeric    | 0                                          | Start position of the selected color palette.                                                                            |
| palette_end       | numeric    | 1                                          | End position of the selected color palette.                                                                              |
| ident_alpha       | numeric    | 0.3                                        | Transparency level for the identity overlay points.                                                                      |
| ident_pointsize   | numeric    | 1                                          | Point size for the identity overlay.                                                                                     |
| filter_ident      | vector     | NULL                                       | Optional; a filter to include only specific identities from group.by.                                                    |
| feature_on_top    | logical    | TRUE                                       | If TRUE, the feature layer (binned data) is rendered on top of the identity overlay.                                     |

---

### `ImageBinPlotFOVs`
Generates bin plots for multiple fields of view (FOVs) in a Seurat object.  
- Optionally enforces the same scale across FOVs using `same_scale`.

| Parameter         | Type       | Default Value                              | Description                                                                                                                |
|-------------------|------------|--------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| obj               | Seurat     | —                                          | A Seurat object containing spatial transcriptomics data.                                                                |
| feature           | character  | —                                          | The gene or feature to be visualized.                                                                                    |
| group.by          | character  | NULL                                       | Optional; column for identity overlay.                                                                                   |
| fovs              | vector     | NULL (defaults to all FOVs via Images(obj))  | A vector of FOV names to plot. If NULL, plots all available FOVs in the object.                                           |
| bin_size          | numeric    | 10                                         | Size of spatial bins (in micrometers).                                                                                   |
| min               | numeric    | 0                                          | Minimum value for the color scale.                                                                                       |
| max_quantile      | character  | "q75"                                      | Quantile string to determine the maximum value for the color scale.                                                      |
| x_min             | numeric    | NULL                                       | Optional minimum x-coordinate; derived from data if not provided.                                                        |
| x_max             | numeric    | NULL                                       | Optional maximum x-coordinate; derived from data if not provided.                                                        |
| y_min             | numeric    | NULL                                       | Optional minimum y-coordinate; derived from data if not provided.                                                        |
| y_max             | numeric    | NULL                                       | Optional maximum y-coordinate; derived from data if not provided.                                                        |
| x_col             | character  | "x"                                        | Column name for x-coordinate values.                                                                                     |
| y_col             | character  | "y"                                        | Column name for y-coordinate values.                                                                                     |
| type              | character  | "expression"                               | Type of data to plot; either "expression" or "molecule".                                                                 |
| assay             | character  | "Xenium"                                   | Assay to use for data extraction.                                                                                        |
| layer             | character  | "count"                                    | Data layer to use for extraction.                                                                                        |
| ident_alpha       | numeric    | 0.3                                        | Transparency level for identity overlay points.                                                                        |
| ident_pointsize   | numeric    | 1                                          | Size of identity overlay points.                                                                                         |
| palette           | character  | c('viridis','inferno','magma','plasma','cividis','mako','rocket','turbo') | Color palette vector for plotting.                           |
| palette_begin     | numeric    | 0                                          | Start point for the color palette.                                                                                       |
| palette_end       | numeric    | 1                                          | End point for the color palette.                                                                                         |
| filter_ident      | vector     | NULL                                       | Optional filter to include only selected identities.                                                                   |
| feature_on_top    | logical    | TRUE                                       | Whether the binned feature layer is drawn on top of the identity layer.                                                  |
| same_scale        | logical    | FALSE                                      | If TRUE, forces all FOV plots to use the same color scale.                                                               |
| ncol              | numeric    | NULL                                       | Optional; the number of columns when arranging the multiple FOV plots in a grid.                                         |

---

### ImageBinPlotObjects
### `ImageBinPlotObjects`
Creates spatial binning plots for multiple Seurat objects, optionally across multiple FOVs, and allows same-scale visualization across objects.
- Optionally enforces the same scale across samples and FOVs using `same_scale`.

| Parameter         | Type       | Default Value                              | Description                                                                                                                |
|-------------------|------------|--------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| obj_list          | list       | —                                          | A named list of Seurat objects to be plotted.                                                                            |
| feature           | character  | —                                          | The gene or feature to visualize across objects.                                                                         |
| group.by          | character  | NULL                                       | Optional; grouping variable for identity overlay across objects.                                                         |
| fovs_list         | list       | NULL (defaults to all FOVs for each object)  | List of FOV vectors corresponding to each object; if NULL, defaults to using all FOVs from each object via Images(.x).        |
| bin_size          | numeric    | 10                                         | Size of the spatial bins (in micrometers).                                                                               |
| min               | numeric    | 0                                          | Minimum value for the color scale.                                                                                       |
| max_quantile      | character  | "q75"                                      | Quantile string to determine the maximum value for the color scale.                                                      |
| x_min             | numeric    | NULL                                       | Optional minimum x-coordinate limit; if NULL will be determined per object/FOV.                                          |
| x_max             | numeric    | NULL                                       | Optional maximum x-coordinate limit; if NULL will be determined per object/FOV.                                          |
| y_min             | numeric    | NULL                                       | Optional minimum y-coordinate limit; if NULL will be determined per object/FOV.                                          |
| y_max             | numeric    | NULL                                       | Optional maximum y-coordinate limit; if NULL will be determined per object/FOV.                                          |
| x_col             | character  | "x"                                        | Column name for x-coordinate values.                                                                                     |
| y_col             | character  | "y"                                        | Column name for y-coordinate values.                                                                                     |
| type              | character  | "expression"                               | Data type to visualize; "expression" or "molecule".                                                                      |
| assay             | character  | "Xenium"                                   | Assay used for data extraction in each Seurat object.                                                                    |
| layer             | character  | "count"                                    | Data layer within the assay to use.                                                                                      |
| ident_alpha       | numeric    | 0.3                                        | Transparency level for identity overlay points.                                                                        |
| ident_pointsize   | numeric    | 1                                          | Size of the identity overlay points.                                                                                     |
| palette           | character  | c('viridis','inferno','magma','plasma','cividis','mako','rocket','turbo') | Color palette vector used for plotting.                         |
| palette_begin     | numeric    | 0                                          | Start point for the color palette.                                                                                       |
| palette_end       | numeric    | 1                                          | End point for the color palette.                                                                                         |
| filter_ident      | vector     | NULL                                       | Optional filter to restrict which identities to plot.                                                                  |
| feature_on_top    | logical    | TRUE                                       | Determines if the binned feature layer is drawn over the identity overlay.                                               |
| same_scale        | logical    | FALSE                                      | If TRUE, applies a uniform color scale across all objects and their FOVs.                                                |
| ncol              | numeric    | NULL                                       | Optional; number of columns when arranging the combined plots in a grid.                                                 |


## Contributing
Submit issues or pull requests to improve functionality or add features.

## License
MIT License - see LICENSE for details.