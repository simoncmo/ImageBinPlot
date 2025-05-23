% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImageBinPlot.R
\name{ImageBinPlotObjects}
\alias{ImageBinPlotObjects}
\title{ImageBinPlotObjects}
\usage{
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
  palette = c("viridis", "inferno", "magma", "plasma", "cividis", "mako", "rocket",
    "turbo"),
  palette_begin = 0,
  palette_end = 1,
  filter_ident = NULL,
  feature_on_top = TRUE,
  same_scale = FALSE,
  ncol = NULL
)
}
\arguments{
\item{obj_list}{A named list of spatial transcriptomics objects to be plotted.}

\item{feature}{The feature to be visualized (e.g., gene expression).}

\item{group.by}{Optional grouping variable for coloring points.}

\item{fovs_list}{A list of FOVs for each object. Defaults to NULL, which uses all FOVs in each object.}

\item{bin_size}{The size of the bins for aggregating data. Default is 10.}

\item{min}{Minimum value for the color scale. Default is 0.}

\item{max_quantile}{Maximum quantile for the color scale. Default is "q75".}

\item{x_min, }{x_max, y_min, y_max Optional limits for the x and y axes.}

\item{x_col, }{y_col Column names for x and y coordinates. Default are "x" and "y".}

\item{type}{The type of data to plot (e.g., "expression"). Default is "expression".}

\item{assay}{The assay to use for data extraction. Default is "Xenium".}

\item{layer}{The data layer to use. Default is "count".}

\item{ident_alpha}{Transparency level for points. Default is 0.3.}

\item{ident_pointsize}{Size of points in the plot. Default is 1.}

\item{palette}{A vector of color palettes to choose from. Default includes 'viridis', 'inferno', etc.}

\item{palette_begin, }{palette_end Range of the color palette to use. Default is 0 to 1.}

\item{filter_ident}{Optional filter for specific identities.}

\item{feature_on_top}{Logical, whether to plot the feature layer on top. Default is TRUE.}

\item{same_scale}{Logical, whether to use the same scale for all objects. Default is FALSE.}

\item{ncol}{Number of columns for arranging plots. Default is NULL.}
}
\value{
A combined plot of all objects and their respective FOVs, optionally with the same scale.
}
\description{
This function generates bin plots for multiple spatial transcriptomics objects, optionally across multiple FOVs.
It allows customization of bin size, color palette, and scaling options. If no FOVs are specified, it defaults to
using all FOVs available in each object.
}
\examples{
# Generate plots for all objects and their FOVs
ImageBinPlotObjects(obj_list = list(obj1 = my_object1, obj2 = my_object2), feature = "GeneA")

# Generate plots for specific FOVs in each object
ImageBinPlotObjects(obj_list = list(obj1 = my_object1, obj2 = my_object2), 
                    feature = "GeneA", 
                    fovs_list = list(c("FOV1", "FOV2"), c("FOV3", "FOV4")))

# Generate plots with the same scale across all objects
ImageBinPlotObjects(obj_list = list(obj1 = my_object1, obj2 = my_object2), 
                    feature = "GeneA", 
                    same_scale = TRUE)

}
\seealso{
\code{\link{ImageBinPlotFOVs}}
}
