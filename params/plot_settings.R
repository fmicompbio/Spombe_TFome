## Colors
## general color map (IBM colors)
main_colors <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
na_color <- "gray70"
bg_color <- "gray90"
pt_color <- "black"
## Binary heatmaps
binary_heatmap_colors <- c("FALSE" = "white", "TRUE" = "darkblue")
enrichment_heatmap_colors <- rev(hcl.colors(9, "RdBu"))
## Complexes
complex_colors <- c(SAGA = "#aa30aa", NuA4 = "#46c8f5")
## Peak set colors
peakset_colors <- c("Common (ubiquitous)" = main_colors[1],
                    "Common peaks (ubiquitous)" = main_colors[1],
                    "Common (frequent)" = main_colors[3],
                    "Common peaks (frequent)" = main_colors[3],
                    "Specific" = main_colors[5],
                    "Specific peaks" = main_colors[5])
