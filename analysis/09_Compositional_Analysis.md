---
title: "Looking at Sample Taxonomic Composition" 
author: "Augustus Pendleton"
date: "12 October, 2026"
knitr: 
  opts_chunk: 
    fig.path: ../figures/09_Compositional_Analysis/
editor_options: 
  chunk_output_type: console
---



# Goals of this document

I am trying to keep fewer steps within each document. So the goal of this document is JUST to use beta diversity measures to produce ordinations and hierarchical dendrograms. There will be a separate document (Compositional Analysis) where I spend time parsing through taxonomic composition of our samples. 

# Load Data and Packages


::: {.cell layout-align="center"}

```{.r .cell-code}
pacman::p_load(phyloseq, ggdendro, dendextend, vegan, patchwork, fun.gus, GUniFrac, tidyverse)

source("code/R/plotting_aesthetics.R")

load("data/08_Beta_Diversity/trim_physeq.RData")

load("data/08_Beta_Diversity/bray_dendro.RData")

trim_melted <- 
  trim_physeq %>%
  psmelt()
```
:::




Now using our grouping by cyanos specifically


::: {.cell layout-align="center"}

```{.r .cell-code}
load("data/08_Beta_Diversity/cyano_physeq.RData")

load("data/08_Beta_Diversity/cy_dendro.RData")

cy_order <- ggdendro::dendro_data(cy_dendro)$label$label

cy_melt <- cyano_physeq %>%
  psmelt()

cy_annotation <- 
  cy_melt %>%
  select(DNA_ID, Cy_Group, Month, Basin) %>%
  unique() %>%
  mutate(DNA_ID = factor(DNA_ID, levels = rev(cy_order)),
         Cy_Group_Name = case_when(
           Cy_Group == 1 ~ "Deep Desert",
           Cy_Group == 2 ~ "Peak Bloom",
           Cy_Group == 3 ~ "Eastern Basin",
           Cy_Group == 4 ~ "Sep. Central Basin",
           Cy_Group == 5 ~ "Bloom Adjacent 2",
           Cy_Group == 6 ~ "Aug. Central Basin",
           Cy_Group == 7 ~ "Bloom Adjacent 1"),
         Cy_Group_Name = factor(Cy_Group_Name,
                                levels = c(
                                  "Deep Desert",
                                  "Peak Bloom",
                                  "Bloom Adjacent 1",
                                  "Bloom Adjacent 2",
                                  "Sep. Central Basin",
                                  "Aug. Central Basin",
                                  "Eastern Basin"
                                ))) %>%
  dplyr::filter(DNA_ID != "KH_D076")

group_plot <- cy_annotation %>%
  ggplot(aes(x=1, y = DNA_ID)) + 
  geom_tile(aes(fill = Cy_Group_Name)) + 
  scale_fill_manual(values = cy_group_colors, name = "Cyano Group") + 
  theme_void()+ 
  theme(legend.position = "left")


cy_bars <- cy_melt %>%
  filter(DNA_ID != "KH_D076") %>%
  mutate(DNA_ID = factor(DNA_ID, levels = rev(cy_order))) %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Family) %>%
  summarize(Fam_Abund = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Family) %>% 
  mutate(Family = ifelse(max(Fam_Abund) > 100000, Family, "(Rare/Unknown)"),
         Family = str_remove(Family, "f__"),
         Family = ifelse(is.na(Family), "(Rare/Unknown)",Family),
         Family = ifelse(Family == "Microcystaceae_A_24132", "Microcystaceae", Family)) %>%
  ungroup() %>% 
  ggplot(aes(x = Fam_Abund, y = DNA_ID, fill = Family)) + 
  geom_col() +
  scale_fill_manual(values = cyano_family_colors) + 
  labs(x = "Absolute Abundance") + 
  theme(axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y= element_blank(),
        plot.margin = margin()) + 
  scale_x_continuous(expand = c(0,0),
                     labels= scales::label_comma()) 

group_plot + cy_bars +
  plot_layout(nrow = 1, 
              widths = c(1,25),
              guides = "collect") & 
  theme(legend.position = "left")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/Figure3b-1.png){fig-align='center' width=2700}
:::

```{.r .cell-code}
# ggsave(filename = "test.png", width = 6, height = 4, dpi = 300)
cy_rel <- cy_melt %>%
  mutate(DNA_ID = factor(DNA_ID, levels = cy_order)) %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Family) %>%
  summarize(Fam_Abund = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Family) %>% 
  mutate(Family = ifelse(max(Fam_Abund) > 5000, Family, "(Rare)"),
         Family = str_remove(Family, "f__")) %>%
  ungroup() %>% 
  group_by(DNA_ID) %>%
  mutate(Fam_Abund = Fam_Abund / sum(Fam_Abund)) %>%
  ggplot(aes(x = Fam_Abund, y = DNA_ID, fill = Family)) + 
  geom_col() +
  scale_fill_manual(values = cyano_family_colors) + 
  labs(x = "Relative Abundance") + 
  theme(axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y= element_blank(),
        plot.margin = margin()) + 
  scale_x_continuous(expand = c(0,0))

group_plot + cy_rel +
  plot_layout(nrow = 1, 
              widths = c(1,16),
              guides = "collect") & 
  theme(legend.position = "left")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/Figure3b-2.png){fig-align='center' width=2700}
:::

```{.r .cell-code}
plot(cy_dendro)
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/Figure3b-3.png){fig-align='center' width=2700}
:::

```{.r .cell-code}
 cy_genus_abs <- 
   cy_melt %>%
  filter(DNA_ID != "KH_D076") %>%
  mutate(DNA_ID = factor(DNA_ID, levels = rev(cy_order))) %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Family, Genus) %>%
  summarize(Fam_Abund = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Genus) %>%
  mutate(Genus = case_when(is.na(Genus) & Family == "f__Cyanobiaceae" ~"Unassigned Cyanobiaceae",
                           is.na(Genus) ~ "Unassigned",
                           max(Fam_Abund) > 50000~ Genus,
                           TRUE ~ "(Rare)"),
         Genus = str_remove(Genus, "g__")
         ) %>% 
  ungroup() %>%
  group_by(Family) %>% 
  mutate(Family = str_remove(Family, "f__"),
         Genus = factor(Genus, levels = names(cyano_genus_colors))) %>%
   ungroup() %>% 
  ggplot(aes(x = Fam_Abund, y = DNA_ID, fill = Genus)) + 
  geom_col() +
  scale_fill_manual(values = cyano_genus_colors) + 
  scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6), expand = c(0,0)) + 
  labs(x = "Absolute Abundance") + 
  theme(axis.title.y = element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y= element_blank(),
        plot.margin = margin()) 
 
  
 group_plot + cy_genus_abs +
  plot_layout(nrow = 1, 
              widths = c(1,16),
              guides = "collect") & 
  theme(legend.position = "left")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/Figure3b-4.png){fig-align='center' width=2700}
:::

```{.r .cell-code}
 tip_annotation_df <- 
  cyano_physeq %>%
  sample_data() %>%
  data.frame() %>%
  select(Month, DNA_ID)

tip_groups <- data.frame(DNA_ID = dendro_data(cy_dendro)$labels$label) %>%
   left_join(tip_annotation_df)

dend_data <- dendro_data(cy_dendro)

cy_dendro <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_point(data = dend_data$labels, aes(x = x, y = y, color = tip_groups$Month), size = 2) + 
  theme_void() +
  coord_flip() + 
  scale_y_reverse() + 
  scale_x_reverse()+
  scale_color_manual(values = month_colors) +
  scale_alpha_manual(values = c(0.8,0,0.8,0.8)) + 
  labs(color = "Month") + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

cy_dendro +  group_plot + cy_genus_abs +
  plot_layout(nrow = 1, 
              widths = c(8,1, 16),
              guides = "collect") & 
  theme(legend.position = "left")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/Figure3b-5.png){fig-align='center' width=2700}
:::

```{.r .cell-code}
ggsave("test.png", width = 9, height= 8)
```
:::




# Mapping Cyanos


::: {.cell layout-align="center"}

```{.r .cell-code}
pacman::p_load(sf, terra, tmap, gstat, install = FALSE)

raw_erie <- read_sf("data/spatial_imports/hydro_p_LakeErie/hydro_p_LakeErie.shp") %>%
  st_transform(crs = "EPSG:4326")

erie_outline <- raw_erie %>%
  st_union()

trim_melted %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID) %>%
  summarize(total_cyanos = sum(Abundance),
            total_phototrophs = mean(total_chl_counts),
            diff = total_cyanos - total_phototrophs) %>%
  pull(diff) %>%
  median()
```

::: {.cell-output .cell-output-stdout}
```
[1] 181693.6
```
:::

```{.r .cell-code}
trim_melted %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID) %>%
  summarize(total_cyanos = sum(Abundance),
            total_phototrophs = mean(total_chl_counts)) %>%
  ggplot(aes(x = total_cyanos, y = total_phototrophs)) + 
  geom_point() + 
  ggpubr::stat_cor(method = "spearman") + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_log10() + 
  scale_x_log10()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/just-cyanos-mapping-1.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
# Cell count??
cy_genus_sf <- 
  trim_melted %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Longitude, Latitude, Month, Phylum, total_chl_counts, Family, Genus) %>%
  summarize(Genus_Abund = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Family = str_remove(Family, "f__"),
         Genus = str_remove(Genus, "g__")) %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326")

cy_fam_sf <- 
  cy_melt %>%
  mutate(DNA_ID = factor(DNA_ID, levels = cy_order)) %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Longitude, Latitude, Month, Family) %>%
  summarize(Fam_Abund = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Family) %>% 
  mutate(Family = ifelse(max(Fam_Abund) > 5000, Family, "(Rare)"),
         Family = str_remove(Family, "f__")) %>%
  ungroup() %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326")

cy_genus_sf %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  group_by(DNA_ID, Month) %>%
  summarize(Abundance = sum(Genus_Abund)) %>%
  ggplot() + 
  geom_sf(aes(color =Abundance/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA) + 
  scale_color_viridis_c(option = "H", name = "All Cyanobacteria\n(million cells/ml)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/just-cyanos-mapping-2.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
cy_genus_sf %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  group_by(DNA_ID, Month) %>%
  summarize(Abundance = mean(total_chl_counts)) %>%
  ggplot() + 
  geom_sf(aes(color =Abundance/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA) + 
  scale_color_viridis_c(option = "H", name = "Cells/ml(millions)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/just-cyanos-mapping-3.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
cy_genus_sf %>%
  dplyr::filter(Family  == "Cyanobiaceae") %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  group_by(DNA_ID, Month) %>%
  summarize(Abundance = sum(Genus_Abund)) %>%
  ggplot() + 
  geom_sf(aes(color =Abundance/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA) + 
  scale_color_viridis_c(option = "H", name = "Cyanobiaceae\n(million cells/ml)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/just-cyanos-mapping-4.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
cy_genus_sf %>%
  dplyr::filter(Genus  == "Microcystis") %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color =Genus_Abund/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA)+
    scale_color_viridis_c(option = "H", name = "Microcystis\n(million cells/ml)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/just-cyanos-mapping-5.png){fig-align='center' width=2100}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
cy_genus_sf %>%
  dplyr::filter(Genus  == "Planktothrix",
                Month == "September") %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color =Genus_Abund/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA)+
      scale_color_viridis_c(option = "H", name = "Cells/ml\n(millions)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/sept-north-shore-1.png){fig-align='center' width=1050}
:::

```{.r .cell-code}
cy_genus_sf %>%
  dplyr::filter(Genus  == "Dolichospermum",
                Month == "September") %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color =Genus_Abund/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA)+
        scale_color_viridis_c(option = "H", name = "Cells/ml\n(millions)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/sept-north-shore-2.png){fig-align='center' width=1050}
:::

```{.r .cell-code}
cy_genus_sf %>%
  dplyr::filter(Genus  == "Pseudanabaena",
              Month == "September") %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color =Genus_Abund/1e6), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA)+ 
          scale_color_viridis_c(option = "H", name = "Cells/ml\n(millions)") + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks= element_blank())
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/sept-north-shore-3.png){fig-align='center' width=1050}
:::

```{.r .cell-code}
cy_fam_sf %>%
  dplyr::filter(Family  == "Cyanobiaceae") %>%
  mutate(Month = factor(Month, levels = c("August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color = Fam_Abund), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA)
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/sept-north-shore-4.png){fig-align='center' width=1050}
:::

```{.r .cell-code}
cy_fam_sf %>%
  dplyr::filter(Family  == "Microcystaceae_A_24132") %>%
  mutate(Month = factor(Month, levels = c("August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color = Fam_Abund), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA)
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/sept-north-shore-5.png){fig-align='center' width=1050}
:::

```{.r .cell-code}
cy_fam_sf %>%
  dplyr::filter(Family  == "Pseudanabaenaceae") %>%
  mutate(Month = factor(Month, levels = c("August","September"))) %>%
  ggplot() + 
  geom_sf(aes(color = Fam_Abund), size = 3, alpha = 1) + 
  facet_wrap(~Month, ncol = 1)+
  geom_sf(data = erie_outline, fill = NA) + 
  scale_color_gradient(low = "grey90", high = "red")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/sept-north-shore-6.png){fig-align='center' width=1050}
:::
:::



# Does higher cyanobacterial counts reduce richness?


::: {.cell layout-align="center"}

```{.r .cell-code}
cy_melt %>%
    mutate(Cy_Group_Name = case_when(
           Cy_Group == 1 ~ "Stn 15 Aug.",
           Cy_Group == 2 ~ "Peak Bloom",
           Cy_Group == 3 ~ "Aug. + Sept. Open",
           Cy_Group == 4 ~ "Sept. Central",
           Cy_Group == 5 ~ "West Basin Bloom Adjacent",
           Cy_Group == 6 ~ "Long Point Desert"
         )) %>%
  ggplot(aes(y = Cy_Group_Name, x = Richness)) + 
  geom_boxplot()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-1.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
cy_melt %>%
    mutate(Cy_Group_Name = case_when(
           Cy_Group == 1 ~ "Stn 15 Aug.",
           Cy_Group == 2 ~ "Peak Bloom",
           Cy_Group == 3 ~ "Aug. + Sept. Open",
           Cy_Group == 4 ~ "Sept. Central",
           Cy_Group == 5 ~ "West Basin Bloom Adjacent",
           Cy_Group == 6 ~ "Long Point Desert"
         )) %>%
  ggplot(aes(y = Cy_Group_Name, x = Shannon)) + 
  geom_boxplot()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-2.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
cy_melt %>%
    mutate(Cy_Group_Name = case_when(
           Cy_Group == 1 ~ "Stn 15 Aug.",
           Cy_Group == 2 ~ "Peak Bloom",
           Cy_Group == 3 ~ "Aug. + Sept. Open",
           Cy_Group == 4 ~ "Sept. Central",
           Cy_Group == 5 ~ "West Basin Bloom Adjacent",
           Cy_Group == 6 ~ "Long Point Desert"
         )) %>%
  ggplot(aes(y = Cy_Group_Name, x = Simpson)) + 
  geom_boxplot()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-3.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
sample_data(trim_physeq) %>%
  data.frame() %>%
  ggplot(aes(x = total_chl_counts / avg_cells_per_ml, y = Richness)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-4.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
sample_data(trim_physeq) %>%
  data.frame() %>%
  ggplot(aes(x = total_chl_counts, y = Richness)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-5.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
sample_data(trim_physeq) %>%
  data.frame() %>%
  ggplot(aes(x = total_chl_counts, y = Simpson)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-6.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
trim_melted %>%
  mutate(is_cy = ifelse(Phylum == "p__Cyanobacteriota", "cy", "nc")) %>%
  group_by(DNA_ID, is_cy, Month) %>%
  summarize(Abund = sum(Abundance)) %>%
  dplyr::filter(!is.na(is_cy)) %>%
  pivot_wider(names_from = is_cy, values_from = Abund) %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot(aes(x = cy/1e6, y = nc/1e6)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  labs(x = "Cyano Abundance (million cells/ml)", y = "Non-Cyano Abundance\n(million cells/ml") + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-7.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
rel_nc_df <- trim_melted %>%
  group_by(DNA_ID) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  dplyr::filter(Phylum != "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Month) %>%
  summarize(Rel_NC = sum(Abundance))

trim_melted %>%
  dplyr::filter(Phylum == "p__Cyanobacteriota") %>%
  group_by(DNA_ID, Month) %>%
  summarize(Abs_C = sum(Abundance)) %>%
  inner_join(rel_nc_df) %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot(aes(x = Abs_C/1e6, y = Rel_NC)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  labs(x = "Absolute Abundance Cyanos (million cells /ml)", y = "% Non-Cyanos") + 
  ggpubr::stat_cor(method = "spearman") 
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-8.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
ggsave(filename = "test.png", width = 9, height = 3)
trim_melted %>%
  group_by(DNA_ID, Month) %>%
  summarize(total_chl_counts = mean(total_chl_counts)) %>%
  inner_join(rel_nc_df) %>%
  mutate(Month = factor(Month, levels = c("May", "August","September"))) %>%
  ggplot(aes(x = total_chl_counts, y = Rel_NC)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  labs(x = "Absolute Abundance Cyanos", y = "Relative Abundance Non-Cyanos") + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-9.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
trim_melted %>%
  mutate(is_cy = ifelse(Phylum == "p__Cyanobacteriota", "cy", "nc")) %>%
  group_by(DNA_ID, is_cy, Month) %>%
  summarize(Richness = sum(Abundance > 0)) %>%
  dplyr::filter(!is.na(is_cy)) %>%
  pivot_wider(names_from = is_cy, values_from = Richness) %>%
  ggplot(aes(x = cy, y = nc)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-10.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
trim_physeq %>%
  sample_data() %>%
  data.frame() %>% 
  #dplyr::filter(Microcystin_Total > 0) %>%
  ggplot(aes(x = Microcystin_Total, y = Richness)) + 
  geom_point()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-11.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
trim_physeq %>%
  sample_data() %>%
  data.frame() %>% 
  dplyr::filter(Microcystin_Total > 0) %>%
  ggplot(aes(x = Microcystin_Total, y = Richness)) + 
  geom_point()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-12.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
trim_melted %>%
  dplyr::filter(Phylum != "p__Cyanobacteriota") %>%
  group_by(DNA_ID) %>%
  summarize(Abs_NC = sum(Abundance)) %>%
  left_join(data.frame(sample_data(trim_physeq))) %>%
  #dplyr::filter(Microcystin_Total > 0) %>%
  ggplot(aes(x = Microcystin_Total, y = Abs_NC)) + 
  geom_point() + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-13.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
trim_melted %>%
  dplyr::filter(Phylum != "p__Cyanobacteriota") %>%
  group_by(DNA_ID) %>%
  summarize(Abs_NC = sum(Abundance)) %>%
  left_join(data.frame(sample_data(trim_physeq))) %>%
  dplyr::filter(Microcystin_Total > 0) %>%
  ggplot(aes(x = Microcystin_Total, y = Abs_NC)) + 
  geom_point() + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-1-14.png){fig-align='center' width=2100}
:::
:::


# DIRTY EUKARYOTIC RICHNESS


::: {.cell layout-align="center"}

```{.r .cell-code}
load("data/04a_Eukaryotic_TreeBuilding/tree_euk_physeq.RData")

rared_euk <- tree_euk_physeq %>%
  prune_samples(sample_sums(.) > 500, .)

non_rare_richness <- rared_euk %>%
  estimate_richness() %>%
  tibble::rownames_to_column("DNA_ID")

rare_richness <- rared_euk %>%
  rarefy_even_depth(500) %>%
  estimate_richness() %>%
  tibble::rownames_to_column("DNA_ID")

inner_join(non_rare_richness, rare_richness, by = "DNA_ID") %>%
  ggplot(aes(x = Observed.x, y = Observed.y)) + 
  geom_point()
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-2-1.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
euk_rich <- 
  rare_richness %>%
  select(DNA_ID, Euk_Richness = Observed)

sample_data(trim_physeq) %>%
  data.frame() %>%
  left_join(euk_rich) %>%
  ggplot(aes(x = total_chl_counts, y = Euk_Richness)) + 
  geom_point() + 
  facet_wrap(~Month, scales = "free") + 
  ggpubr::stat_cor(method = "spearman")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-2-2.png){fig-align='center' width=2100}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
cy_melt %>%
  dplyr::filter(DNA_ID %in% c("KH_D064", "KH_D032"),
                Genus == "g__Microcystis") %>%
  ggplot(aes(x = Abundance, y = DNA_ID, fill = ASVseqs)) + 
  geom_col() + 
  theme(legend.position = "none")
```

::: {.cell-output-display}
![](../figures/09_Compositional_Analysis/unnamed-chunk-3-1.png){fig-align='center' width=2100}
:::

```{.r .cell-code}
tax_table(tree_meta_physeq) %>%
  data.frame() %>%
  select(Supergroup:Species) %>%
  unique() %>% 
  arrange(Class) %>%
  write_csv("Metazoa_Detections.csv")
```

::: {.cell-output .cell-output-error}
```
Error in `h()`:
! error in evaluating the argument 'x' in selecting a method for function 'unique': error in evaluating the argument 'object' in selecting a method for function 'tax_table': object 'tree_meta_physeq' not found
```
:::
:::
