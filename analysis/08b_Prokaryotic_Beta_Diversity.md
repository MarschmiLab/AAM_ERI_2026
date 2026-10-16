---
title: "Comparing Samples with Beta Diversity" 
author: "Augustus Pendleton"
date: "12 October, 2026"
knitr: 
  opts_chunk: 
    fig.path: ../figures/08b_Prokaryotic_Beta_Diversity/
editor_options: 
  chunk_output_type: console
---



# Goals of this document

Use dissimilarity metrics (beta diversity) to find patterns or groups in the bacterial composition between our samples. While this seems kinda simple going in, there are a lot of choices to make with incorporating different dissimilarity metrics and factoring in absolute abundance. 

# Load Data and Packages


::: {.cell layout-align="center"}

```{.r .cell-code}
pacman::p_load(phyloseq, ggdendro, dendextend, vegan, patchwork, fun.gus, GUniFrac, tidyverse, install = FALSE)


# Make sure we're using the MarschmiLab version, which incorporates absolute abundance!
packageDescription("GUniFrac")
```

::: {.cell-output .cell-output-stdout}
```
Package: GUniFrac
Type: Package
Title: Generalized UniFrac Distances, Distance-Based Multivariate
        Methods and Feature-Based Univariate Methods for Microbiome
        Data Analysis
Version: 1.8.1
Date: 2023-09-13
Author: Jun Chen, Xianyang Zhang, Lu Yang, Lujun Zhang
Maintainer: Jun Chen <chen.jun2@mayo.edu>
Description: A suite of methods for powerful and robust microbiome data
        analysis including data normalization, data simulation,
        community-level association testing and differential abundance
        analysis. It implements generalized UniFrac distances,
        Geometric Mean of Pairwise Ratios (GMPR) normalization,
        semiparametric data simulator, distance-based statistical
        methods, and feature-based statistical methods. The
        distance-based statistical methods include three extensions of
        PERMANOVA: (1) PERMANOVA using the Freedman-Lane permutation
        scheme, (2) PERMANOVA omnibus test using multiple matrices, and
        (3) analytical approach to approximating PERMANOVA p-value.
        Feature-based statistical methods include linear model-based
        methods for differential abundance analysis of zero-inflated
        high-dimensional compositional data.
Depends: R (>= 3.5.0)
Suggests: ade4, knitr, markdown, ggpubr
Imports: Rcpp (>= 0.12.13), vegan, ggplot2, matrixStats, Matrix, ape,
        parallel, stats, utils, statmod, rmutil, dirmult, MASS,
        ggrepel, foreach, modeest, inline, methods
LinkingTo: Rcpp
NeedsCompilation: yes
VignetteBuilder: knitr
License: GPL-3
Encoding: UTF-8
Packaged: 2026-10-07 20:33:45 UTC; arp277
Repository: CRAN
Date/Publication: 2023-09-14 00:02:32 UTC
RemoteType: github
RemoteHost: api.github.com
RemoteRepo: GUniFrac
RemoteUsername: MarschmiLab
RemoteRef: HEAD
RemoteSha: 91321d0d07a1584b962638f21584be72c6d15ff9
GithubRepo: GUniFrac
GithubUsername: MarschmiLab
GithubRef: HEAD
GithubSHA1: 91321d0d07a1584b962638f21584be72c6d15ff9
Built: R 4.3.3; x86_64-pc-linux-gnu; 2026-10-07 20:33:46 UTC; unix

-- File: /home/arp277/R/x86_64-pc-linux-gnu-library/4.3/GUniFrac/Meta/package.rds 
```
:::

```{.r .cell-code}
source("code/R/plotting_aesthetics.R")

load("data/07_Alpha_Diversity/diversity_physeq.RData")
```
:::



# Actual stuff I need




::: {.cell layout-align="center"}

```{.r .cell-code}
cell_counts_prep1 <- 
  diversity_physeq %>%
  sample_data %>%
  data.frame %>%
  select(DNA_ID, avg_cells_per_ml)

cell_counts1 <- cell_counts_prep1$avg_cells_per_ml

names(cell_counts1) = cell_counts_prep1$DNA_ID

rel_matrix <- diversity_physeq %>%
    transform_sample_counts(., \(x)x/sum(x)) %>%
    otu_table() %>%
    as.matrix()
  
  colnames(rel_matrix)
```

::: {.cell-output .cell-output-stdout}
```
 [1] "KH_D001" "KH_D009" "KH_D017" "KH_D025" "KH_D032" "KH_D040" "KH_D056"
 [8] "KH_D063" "KH_D071" "KH_D079" "KH_D002" "KH_D010" "KH_D018" "KH_D026"
[15] "KH_D033" "KH_D041" "KH_D057" "KH_D064" "KH_D003" "KH_D011" "KH_D019"
[22] "KH_D027" "KH_D034" "KH_D042" "KH_D058" "KH_D065" "KH_D004" "KH_D012"
[29] "KH_D020" "KH_D028" "KH_D035" "KH_D043" "KH_D051" "KH_D059" "KH_D066"
[36] "KH_D074" "KH_D005" "KH_D013" "KH_D021" "KH_D029" "KH_D036" "KH_D044"
[43] "KH_D052" "KH_D060" "KH_D067" "KH_D075" "KH_D083" "KH_D006" "KH_D022"
[50] "KH_D030" "KH_D037" "KH_D045" "KH_D053" "KH_D061" "KH_D076" "KH_D007"
[57] "KH_D015" "KH_D023" "KH_D031" "KH_D038" "KH_D046" "KH_D054" "KH_D062"
[64] "KH_D069" "KH_D077" "KH_D008" "KH_D016" "KH_D039" "KH_D047" "KH_D055"
[71] "KH_D070" "KH_D078"
```
:::

```{.r .cell-code}
  abs_matrix <- map2(colnames(rel_matrix), cell_counts1, \(x,y)round(rel_matrix[,x]*y)) %>%
    bind_cols() %>%
    as.matrix()
  
  row.names(abs_matrix) <- row.names(rel_matrix)
  
  abs_physeq <- diversity_physeq
  
  abs_physeq@otu_table <- NULL
  
  abs_physeq@otu_table <- otu_table(abs_matrix, taxa_are_rows = TRUE)

samples_to_keep <- 
  sample_data(diversity_physeq) %>%
  data.frame() %>%
  filter(Depth_Class == "S") %>%
  pull(DNA_ID)


trim_physeq <- 
  abs_physeq %>%
  subset_samples(., DNA_ID %in% samples_to_keep)

cyano_physeq <- trim_physeq %>%
  subset_samples(., DNA_ID %in% samples_to_keep) %>%
  subset_samples(., Month != "May") %>%
  subset_samples(., DNA_ID != "KH_D076") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  subset_taxa(., Phylum == "p__Cyanobacteriota")

# cyano_prep_physeq <- 
#   diversity_physeq %>%
#   subset_samples(., DNA_ID %in% samples_to_keep) %>%
#   subset_samples(., Month != "May") %>%
#   subset_samples(., DNA_ID != "KH_D076") %>%
#   prune_taxa(taxa_sums(.) > 0, .)

# cell_counts_prep <- 
#   cyano_prep_physeq %>%
#   sample_data %>%
#   data.frame %>%
#   select(DNA_ID, avg_cells_per_ml)
# 
# cell_counts <- cell_counts_prep$avg_cells_per_ml
# 
# names(cell_counts) = cell_counts_prep$DNA_ID
# 
# cyano_bray_rarefaction <-
#   
#   map(c(1:1000), \(x){
#   
#   rare_physeq <- cyano_prep_physeq %>%
#     rarefy_even_depth(rngseed = x, replace = FALSE)
#   
#  
#   
#   rel_matrix <- rare_physeq %>%
#     transform_sample_counts(., \(x)x/sum(x)) %>%
#     otu_table() %>%
#     as.matrix()
#   
#   
#   colnames(rel_matrix)
#   
#   abs_matrix <- map2(colnames(rel_matrix), cell_counts, \(x,y)round(rel_matrix[,x]*y)) %>%
#     bind_cols() %>%
#     as.matrix()
#   
#   row.names(abs_matrix) <- row.names(rel_matrix)
#   
#   abs_rare_physeq <- rare_physeq
#   
#   abs_rare_physeq@otu_table <- NULL
#   
#   abs_rare_physeq@otu_table <- otu_table(abs_matrix, taxa_are_rows = TRUE)
#   
#   
#   rel_cy_physeq <- 
#     rare_physeq %>%
#     subset_taxa(., Phylum == "p__Cyanobacteriota")
#   
#   abs_cy_physeq <-
#     abs_rare_physeq %>%
#     subset_taxa(., Phylum == "p__Cyanobacteriota")
#   
#   
#   rel_cy_bray <- distance(rel_cy_physeq, method = "bray")
#   
#   abs_cy_bray <- distance(abs_cy_physeq, method = "bray")
#   
#   results <- list(Rel_Cy_Bray = rel_cy_bray,
#                   Abs_Cy_Bray = abs_cy_bray)
#   
#   return(results)
#   
# })
# 
# save(cyano_bray_rarefaction, file = "data/08_Beta_Diversity/cyano_bray_rarefaction.RData")

load("data/08_Beta_Diversity/cyano_bray_rarefaction.RData")

rel_dists <- map(cyano_bray_rarefaction, \(x)x$Rel_Cy_Bray)
abs_dists <- map(cyano_bray_rarefaction, \(x)x$Abs_Cy_Bray)


rel_bray_dist <- Reduce(`+`, rel_dists) / length(rel_dists)
abs_bray_dist <- Reduce(`+`, abs_dists) / length(abs_dists)
```
:::



# Clustering


::: {.cell layout-align="center"}

```{.r .cell-code}
rel_clust <- hclust(rel_bray_dist, method = "average")
abs_clust <- hclust(abs_bray_dist, method = "average")

dend_list <- dendlist(as.dendrogram(rel_clust), as.dendrogram(abs_clust))

dend_list_rotated <- untangle(dend_list, method = "step2side")

# Plot the tanglegram
tanglegram(dend_list_rotated, 
           common_subtrees_color_lines = TRUE, 
           highlight_distinct_edges = TRUE)
```

::: {.cell-output-display}
![](../figures/08b_Prokaryotic_Beta_Diversity/unnamed-chunk-1-1.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
# Measure entanglement (0 = perfect alignment, 1 = high tangle)
entanglement(dend_list_rotated)
```

::: {.cell-output .cell-output-stdout}
```
[1] 0.004503066
```
:::

```{.r .cell-code}
sort(sample_sums(cyano_prep_physeq))
```

::: {.cell-output .cell-output-error}
```
Error in `h()`:
! error in evaluating the argument 'object' in selecting a method for function 'otu_table': object 'cyano_prep_physeq' not found
```
:::
:::



Very similar results using relative versus absolute abundances!



::: {.cell layout-align="center"}

```{.r .cell-code}
library(sf)
library(terra)

raw_erie <- read_sf("data/spatial_imports/hydro_p_LakeErie/hydro_p_LakeErie.shp") %>%
  st_transform(crs = "EPSG:4326")

erie_outline <- raw_erie %>%
  st_union()

rel_cuts <- cutree(rel_clust, 7)  

abs_cuts <- cutree(abs_clust, 7)


cy_sf <- cyano_physeq %>%
  sample_data() %>%
  data.frame() %>%
  left_join(data.frame(DNA_ID = names(rel_cuts), Rel_Cy_Group = factor(rel_cuts))) %>%
  left_join(data.frame(DNA_ID = names(abs_cuts), Abs_Cy_Group = factor(abs_cuts))) %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326")

cy_sf %>%
  ggplot() + 
  geom_sf(data = erie_outline, fill = NA) + 
  geom_sf(aes(color = Rel_Cy_Group), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)
```

::: {.cell-output-display}
![](../figures/08b_Prokaryotic_Beta_Diversity/unnamed-chunk-2-1.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
cy_sf %>%
  ggplot() + 
  geom_sf(data = erie_outline, fill = NA) + 
  geom_sf(aes(color = Abs_Cy_Group), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)
```

::: {.cell-output-display}
![](../figures/08b_Prokaryotic_Beta_Diversity/unnamed-chunk-2-2.png){fig-align='center' width=2100}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
cy_sf %>%
  mutate(Month = factor(Month, levels = c("May", "August","September")),
         Cy_Group_Name = case_when(
           Abs_Cy_Group == 1 ~ "Deep Desert",
           Abs_Cy_Group == 2 ~ "Peak Bloom",
           Abs_Cy_Group == 3 ~ "Eastern Basin",
           Abs_Cy_Group == 4 ~ "Sep. Central Basin",
           Abs_Cy_Group == 5 ~ "Bloom Adjacent 2",
           Abs_Cy_Group == 6 ~ "Aug. Central Basin",
           Abs_Cy_Group == 7 ~ "Bloom Adjacent 1"),
         Cy_Group_Name = factor(Cy_Group_Name,
                                levels = c(
                                  "Peak Bloom",
                                  "Bloom Adjacent 1",
                                  "Bloom Adjacent 2",
                                  "Eastern Basin",
                                  "Aug. Central Basin",
                                  "Sep. Central Basin",
                                  "Deep Desert"
                                ))) %>%
  ggplot() + 
  geom_sf(data = erie_outline, fill = NA) + 
  geom_sf(aes(color = Cy_Group_Name), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  scale_color_manual(values = cy_group_colors, name = "Cyano Group") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank(),
        legend.position = "none")
```

::: {.cell-output-display}
![](../figures/08b_Prokaryotic_Beta_Diversity/Figure3a-1.png){fig-align='center' width=900}
:::

```{.r .cell-code}
ggsave(filename = "test.png", height = 8, width = 8, dpi = 300)




cyano_physeq <- cyano_physeq %>%
  microViz::ps_join(data.frame(DNA_ID = names(abs_cuts), Cy_Group = factor(abs_cuts)))

df <- plot_ordination(cyano_physeq, 
                      ordination = ordinate(cyano_physeq, distance = "bray", method = "PCoA"),
                      color = "Cy_Group", justDF = TRUE) 

df %>%
  mutate( Cy_Group_Name = case_when(
    Cy_Group == 1 ~ "Deep Desert",
    Cy_Group == 2 ~ "Peak Bloom",
    Cy_Group == 3 ~ "Eastern Basin",
    Cy_Group == 4 ~ "Sep. Central Basin",
    Cy_Group == 5 ~ "Bloom Adjacent 2",
    Cy_Group == 6 ~ "Aug. Central Basin",
    Cy_Group == 7 ~ "Bloom Adjacent 1")
  ) %>%
  ggplot(aes(x = Axis.1, y = Axis.2, color = Cy_Group_Name)) + 
  geom_point() + 
  geom_point(size = 4, alpha = 0.4) + 
  scale_color_manual(values = cy_group_colors) + 
  theme(legend.position = "none") + 
  labs(x = "PCoA 1 [25.6%]", y = "PCoA 2 [20.9%]", title = "Relative Unifrac")
```

::: {.cell-output-display}
![](../figures/08b_Prokaryotic_Beta_Diversity/Figure3a-2.png){fig-align='center' width=900}
:::

```{.r .cell-code}
ggsave("test.png", width = 3, height = 3)

save(cyano_physeq, file = "data/08_Beta_Diversity/cyano_physeq.RData")

cy_dendro <- abs_clust %>%
  as.dendrogram()

save(cy_dendro, file = "data/08_Beta_Diversity/cy_dendro.RData")
```
:::




# Saving Objects


::: {.cell layout-align="center"}

```{.r .cell-code}
save(abs_physeq,
     file = "data/08_Beta_Diversity/abs_physeq.RData")

save(trim_physeq,
     file = "data/08_Beta_Diversity/trim_physeq.RData")
```
:::
