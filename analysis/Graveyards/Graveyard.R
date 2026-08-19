
filt_split <- cyano_melted %>%
  left_join(top_cyano_asvs) %>% 
  rowwise() %>%
  filter(ASV %in% ASVs) %>% 
  group_by(Month) %>%
  group_split()

names(filt_split) <- c("August","May","September")

trans_mats <- 
  map(filt_split, \(x){
    mat1 <- 
      x %>%
      select(OTU, Sample, Abundance) %>%
      pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0) %>%
      column_to_rownames("Sample") %>%
      as.matrix()
    
    log_mat <- log10(mat1 + 1)
  })

trans_asv_vals <- 
  map(trans_mats, asplit, 2)

dt_results <- 
  map(trans_asv_vals, \(x){
    dt_res <- map(x, diptest::dip.test, simulate.p.value = TRUE, B=2000)
    ps <- map_dbl(dt_res, \(x)x$p.value)
    bimod_asvs <- names(ps)[ps<0.01]
  })

dt_results

map2(filt_split[c(1,3)], dt_results[c(1,3)], \(x,y){
  x %>%
    filter(ASV %in% y) %>%
    ggplot(aes(x = log10(Abundance), group = ASV)) + 
    geom_density() +
    facet_wrap(~Genus, scales = "free")
})

# Okay, those seem reasonable!


unimod_mats <- 
  map2(trans_mats, dt_results, \(x,y)x[, !(colnames(x) %in% y)])

unimod_asv_vals <- 
  map(unimod_mats, asplit, 2)

dist_names <- c("Normal", "Logistic", "Exponential")
dists_to_test <- c("norm","logis","exp")

names(dists_to_test) <- dist_names

test_fit <- function(name, values, dists){
  aics <- map_dbl(dists, \(dist){
    
    model <- tryCatch(
      expr = {
        fitdistrplus::fitdist(as.numeric(values), dist)
      },
      error = function(e){list(aic = Inf)}
    )
    
    return(model$aic)
  })
  
  as_tibble_row(aics) %>%
    mutate(ASV = name)
}

model_fits <- 
  map(unimod_asv_vals, \(uav){
    map2(uav, names(uav), 
         \(x,y)test_fit(name = y, values = x, dists = dists_to_test)) %>%
      list_rbind(names_to = "ASV")
  }
  )

top_models <- 
  map(model_fits, \(x){
    x %>%
      pivot_longer(Normal:Exponential, names_to = "Model", values_to = "AIC") %>%
      group_by(ASV) %>%
      slice_min(AIC) %>%
      ungroup() %>%
      dplyr::select(-AIC)
  })

map(top_models, \(x)count(x, Model))


bi_df <- 
  map(dt_results[c(1,3)], \(x){
    data.frame(ASV = x,
               Model = "Bimodal")
  }) %>%
  bind_rows(.id = "Month")

all_mod <-
  top_models %>%
  bind_rows(.id = "Month") %>%
  rbind(bi_df) %>%
  left_join(sum_cyanos_asv)


all_mod %>%
  count(Month, Model)

all_mod %>%
  filter(Month == "May")

all_mod %>%
  filter(Month == "August") %>%
  count(Genus, Model) %>%
  ggplot(aes(x = n, y = Genus, fill = Model)) + 
  geom_col()

all_mod %>%
  filter(Month == "September") %>%
  count(Genus, Model) %>%
  ggplot(aes(x = n, y = Genus, fill = Model)) + 
  geom_col()


big_genera <- cyano_melted %>%
  filter(!is.na(Genus)) %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Total_Abund = sum(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Median_Abund = median(Abundance)) %>%
  ungroup() %>%
  group_by(Month) %>%
  slice_max(n = 20, order_by = Total_Abund) %>%
  pull(Genus) %>%
  unique()

cyano_melted %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Total_Abund = sum(Abundance)) %>%
  mutate(Month = factor(Month, levels = c("May","August","September"))) %>%
  ungroup() %>%
  mutate(Genus = ifelse(Genus %in% big_genera,
                        Genus, 
                        "(Rare)"),
         Genus = str_remove(Genus, "g__")) %>%
  ggplot(aes(x = Prevalence, y = Max_Abund, color = Genus, size = Total_Abund)) + 
  geom_point(position = position_jitter(height = 0, width = 0.01),
             alpha = 0.7) + 
  scale_y_log10()+ 
  facet_wrap(~Month) + 
  scale_color_manual(values = cyano_genus_colors)

cyano_melted %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Total_Abund = sum(Abundance)) %>%
  mutate(Month = factor(Month, levels = c("May","August","September"))) %>%
  ungroup() %>%
  mutate(Genus = ifelse(Genus %in% big_genera,
                        Genus, 
                        "(Rare)"),
         Genus = str_remove(Genus, "g__")) %>%
  filter(Genus != "(Rare)") %>%
  ggplot(aes(x = Prevalence, y = Max_Abund, color = Genus, size = Total_Abund)) + 
  geom_point(position = position_jitter(height = 0, width = 0.01),
             alpha = 0.7) + 
  scale_y_log10()+ 
  facet_wrap(~Month) + 
  scale_color_manual(values = cyano_genus_colors)


centroids <- cyano_melted %>%
  mutate(Genus = ifelse(Genus %in% big_genera,
                        Genus, 
                        "(Rare)"),
         Genus = str_remove(Genus, "g__")) %>%
  filter(Genus != "(Rare)") %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Total_Abund = sum(Abundance)) %>%
  mutate(Month = factor(Month, levels = c("May","August","September"))) %>%
  ungroup() %>%
  group_by(Month, Genus) %>%
  summarize(y_centr = sum(Max_Abund * Total_Abund) / sum(Total_Abund),
            x_centr = sum(Prevalence * Total_Abund) / sum(Total_Abund),
            Total_Abund = sum(Total_Abund)) %>%
  ungroup() %>%
  rename(Prevalence = x_centr, Max_Abund = y_centr)




cyano_melted %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Total_Abund = sum(Abundance)) %>%
  mutate(Month = factor(Month, levels = c("May","August","September"))) %>%
  ungroup() %>%
  mutate(Genus = ifelse(Genus %in% big_genera,
                        Genus, 
                        "(Rare)"),
         Genus = str_remove(Genus, "g__")) %>%
  filter(Genus != "(Rare)") %>%
  ggplot(aes(x = Prevalence, y = Max_Abund, color = Genus, size = Total_Abund)) + 
  geom_point(position = position_jitter(height = 0, width = 0.01),
             alpha = 0.2) + 
  geom_point(data = centroids, size = 5) + 
  scale_y_log10()+ 
  facet_wrap(~Month) + 
  scale_color_manual(values = cyano_genus_colors)


cyano_melted %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Total_Abund = sum(Abundance)) %>%
  mutate(Month = factor(Month, levels = c("May","August","September"))) %>%
  ungroup() %>%
  mutate(Genus = ifelse(Genus %in% big_genera,
                        Genus, 
                        "(Rare)"),
         Genus = str_remove(Genus, "g__")) %>%
  filter(Genus != "(Rare)") %>%
  ggplot(aes(x = Prevalence, y = Max_Abund, color = Genus, size = Total_Abund)) + 
  geom_point(position = position_jitter(height = 0, width = 0.01),
             alpha = 0) + 
  geom_point(data = centroids, size = 5) + 
  scale_y_log10()+ 
  facet_wrap(~Month) + 
  scale_color_manual(values = cyano_genus_colors)

centroids %>% 
  ggplot(aes(x = Prevalence, y = Max_Abund)) + 
  geom_path(aes(group = Genus), arrow = arrow()) + 
  #geom_point(aes(size = log10(Total_Abund), shape = Month)) +
  scale_y_log10() + 
  #scale_shape_manual(values = c(17,15,16)) +
  facet_wrap(~Genus)

cyano_melted %>%
  group_by(Month, ASV, Genus, Species) %>%
  summarize(Max_Abund = max(Abundance),
            Prevalence = sum(Abundance>0)/n(),
            Total_Abund = sum(Abundance)) %>%
  mutate(Month = factor(Month, levels = c("May","August","September"))) %>%
  ungroup() %>%
  mutate(Genus = ifelse(Genus %in% big_genera,
                        Genus, 
                        "(Rare)"),
         Genus = str_remove(Genus, "g__")) %>%
  filter(Genus != "(Rare)") %>%
  arrange(Month, Genus) %>%
  ggplot(aes(x = Prevalence, y = Max_Abund)) + 
  geom_path(aes(group = ASV), alpha = 0.2)+
  geom_point(aes(color = Month)) + 
  scale_y_log10() + 
  facet_wrap(~Genus)


library(sf)

geo_dists <- map(filt_split, \(x){
  
  sf_df <- 
    x %>%
    select(Rep_ID, Latitude, Longitude) %>%
    unique() %>%
    st_as_sf(coords = c("Longitude","Latitude"), crs = "EPSG:4326") %>%
    st_transform(crs = "EPSG:3174")
  
  geo_dist_mat <- 
    sf_df %>%
    st_distance()
  
  row.names(geo_dist_mat) <- sf_df$Rep_ID
  colnames(geo_dist_mat) <- sf_df$Rep_ID
  
  diag(geo_dist_mat) <- NA
  
  geo_dist_mat[upper.tri(geo_dist_mat)] <- NA
  
  geo_mat_df <- as.data.frame(geo_dist_mat)
  
  geo_mat_df %>%
    mutate(item1 = rownames(geo_mat_df)) %>%
    pivot_longer(!item1, names_to = "item2", values_to = "geo_dist") %>%
    filter(!is.na(geo_dist))
})



asv_dists <- 
  map(filt_split, \(month_df) {
    grouped <- month_df %>%
      select(ASV, Rep_ID, Abundance) %>%
      group_by(ASV)
    
    split <- group_split(grouped)
    
    names(split) <- group_keys(grouped)$ASV
    
    tidy_dists <- map(split, \(x) {
      dist_mat <-
        x %>%
        select(-ASV) %>%
        column_to_rownames("Rep_ID") %>%
        dist() %>%
        as.matrix()
      
      
      diag(dist_mat) <- NA
      
      dist_mat[upper.tri(dist_mat)] <- NA
      
      dist_df <- as.data.frame(dist_mat)
      
      dist_df %>%
        mutate(item1 = rownames(dist_df)) %>%
        pivot_longer(!item1, names_to = "item2", values_to = "asv_dist") %>%
        filter(!is.na(asv_dist))
      
    })
  })

asv_cor_res <- map2(asv_dists, geo_dists, \(ads, gd){
  
  map(ads, \(ad){
    ad %>%
      left_join(gd) %>%
      rstatix::cor_test(asv_dist, geo_dist, method = "spearman")
  }) %>%
    bind_rows(.id = "ASV")
  
})  %>%
  bind_rows(.id = "Month")





asv_cor_res %>%
  ggplot(aes(x = cor, y=-log10(p), color = Month)) + 
  geom_point() 

asv_cor_res %>%
  left_join(asv_ids) %>%
  mutate(Genus = str_remove(Genus, "g__")) %>%
  ggplot(aes(x = cor, y=-log10(p), color = Genus)) + 
  geom_point() +
  facet_wrap(~Month) + 
  scale_color_manual(values = cyano_genus_colors)


sig_asvs <- asv_cor_res %>%
  filter(p < 0.01, abs(cor) > 0.2) %>%
  select(Month, ASV) %>%
  group_by(Month) %>%
  group_split()

asv_dists[[1]]

sig_asv_dists <- 
  
  map2(asv_dists, sig_asvs, \(x,y){
    
    x[names(x) %in% y$ASV]
    
  })

asv_ids

library(units)

plots <- map2(sig_asv_dists, geo_dists, \(sads, gd){
  
  map(sads, \(ad){
    ad %>%
      left_join(gd) %>%
      ggplot(aes(x = geo_dist, y = -asv_dist)) + 
      geom_point() + 
      geom_smooth(method = "lm", se = FALSE, color = "black") 
    
  })
  
})


names(plots)
map(plots, names)

map(plots, \(x)x)
