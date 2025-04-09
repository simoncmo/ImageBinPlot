# ImageBinPlot

A Simple R function for binning and visualizing spatial transcriptomics data from Seurat objects. This utility supports both cell-based expression data and molecule-based transcript counts, with optional identity overlays using `scattermore`. Designed for Xenium data but adaptable to other spatial datasets.

## Features
- **Binning**: Bin spatial data into a grid with customizable bin sizes.
- **Visualization**: Generate raster plots with `ggplot2` and the `viridis` color scale.
- **Flexibility**: Supports two data types:
  - `expression`: Cell-level transcript counts.
  - `molecule`: Individual transcript coordinates.
- **Identity Overlay**: Overlay cell identities (e.g., clusters) using `geom_scattermore`.
- **Customizable**: Adjust bin size, color scaling, transparency, and layer order.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/simoncmo/ImageBinPlot.git
   ```

2. **Source the Script**:
   In your R session, source the script directly:
   ```R
   source("ImageBinPlot.R")
   ```

## Dependencies

Install the required R packages:
```R
install.packages(c("tidyverse", "Seurat", "ggplot2", "glue", "scattermore"))
```

- `tidyverse`: For data manipulation and visualization (`dplyr`, `ggplot2`, etc.).
- `Seurat`: For handling spatial transcriptomics data.
- `ggplot2`: For plotting (included in `tidyverse`).
- `glue`: For string interpolation in plot labels.
- `scattermore`: For efficient scatter plot overlays.

## Usage

### Basic Setup
Load your Seurat object (e.g., Xenium data) and source the script:
```R
obj <- readRDS("/path/to/your_seurat_object.rds")
library(qs)  # Optional, for loading .qs files
obj <- qread("/path/to/your_seurat_object.qs", nthreads = 30)
source("ImageBinPlot.R")
```

### Example 1: Expression-Based Binning
Plot binned expression data for a gene:
```R
# Plot expression of "EPCAM" with 10µm bins
p <- ImageBinPlot(
  obj, 
  feature = "EPCAM", 
  bin_size = 10, # in µm
  max_multiplier = 0.3, 
  type = "expression"
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
  max_multiplier = 0.15, 
  type = "molecule"
)
ggsave("tp63_molecules.pdf", p, width = 10, height = 10)
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
  max_multiplier = 0.15, 
  type = "expression",
  ident_alpha = 0.3, 
  ident_pointsize = 1, 
  feature_on_top = TRUE
)
ggsave("epcam_with_clusters.pdf", p, width = 10, height = 10)
```

### Example 4: Filtered Identity Overlay
Show only specific clusters:
```R
# Plot "EPCAM" with only cluster "0" overlaid. Can provide more as a vector. e.g. c("0","10")
p <- ImageBinPlot(
  obj, 
  feature = "EPCAM", 
  group.by = "seurat_clusters", 
  filter_ident = "0", 
  bin_size = 10, 
  max_multiplier = 0.3, 
  type = "expression",
  ident_alpha = 0.8, 
  feature_on_top = TRUE
)
ggsave("EPCAM_cluster0.pdf", p, width = 10, height = 10)
```

### Example 5: Batch Processing
Generate plots for multiple genes:
```R
genes <- c("EPCAM", "TP63", "EPCAM")
pdf("batch_plots.pdf", width = 10, height = 10)
for (gene in genes) {
  p <- ImageBinPlot(
    obj, 
    feature = gene, 
    bin_size = 10, 
    max_multiplier = 0.3, 
    type = "expression"
  )
  print(p)
}
dev.off()
```

## Function Details

### `ImageBinPlot`
Main function for binning and plotting spatial data.

**Arguments**:
- `obj`: Seurat object with spatial data.
- `feature`: Gene name to plot.
- `group.by`: Optional identity column (e.g., "seurat_clusters").
- `type`: `"expression"` (cell-based) or `"molecule"` (transcript-based). Default: `"expression"`.
- `bin_size`: Size of bins in µm. Default: 10.
- `max_multiplier`: Scaling factor for color limits. Default: 0.3.
- `min`: Minimum value for color scale. Default: 0.
- `ident_alpha`: Transparency of identity points. Default: 0.3.
- `ident_pointsize`: Size of identity points. Default: 1.
- `filter_ident`: Optional vector of identities to filter. Default: NULL.
- `feature_on_top`: If TRUE, feature raster is above identity points. Default: TRUE.

**Returns**: A `ggplot` object.

## Notes
- Ensure your Seurat object contains the required spatial data (e.g., centroids for expression, molecules for transcripts).
- Adjust `max_multiplier` to fine-tune color scaling based on your data range.
- The script assumes Xenium data structure but can be adapted for other formats.

## Contributing
Feel free to submit issues or pull requests to improve functionality or add features!

## License
MIT License - see [LICENSE](LICENSE) for details.
