toxin_plotter <- function(df_in_nest_raw, 
                            title_lbl, 
                            var_limits, 
                            coop, 
                            raster_grid, 
                            certainty_station_list, 
                            var_of_interest, 
                            legend_title, 
                            palette = "tol.sunset", 
                            legend_show = TRUE,
                            int_deg = 3){
  
  df_in_nest <- df_in_nest_raw[!is.na(df_in_nest_raw[[var_of_interest]]),]
  
  cert_df <- certainty_station_list %>%
    mutate(certainty = ifelse(Station_ID %in% df_in_nest$Station_ID, 1, 0))
  
  some_missing <- any(cert_df$certainty == 0)
  # Loading spatial data
  
  res <- gstat(formula = as.formula(paste0(var_of_interest, "~1")), locations = df_in_nest, 
               nmax = nrow(df_in_nest), 
               set = list(idp = int_deg))
  
  resp <- predict(res, coop) %>% st_sf()
  
  resp$x <- st_coordinates(resp)[,1]
  resp$y <- st_coordinates(resp)[,2]
  resp$pred <- resp$var1.pred
  
  output_raster <- terra::rasterize(resp, raster_grid, field = "pred", fun = "mean")
  
  
  if(some_missing){
    res_cert <- gstat(formula = certainty ~ 1, locations = cert_df, 
                      nmax = nrow(cert_df), 
                      set = list(idp = 5))
    
    resp_cert <- predict(res_cert, coop)
    
    resp_cert$x <- st_coordinates(resp_cert)[,1]
    resp_cert$y <- st_coordinates(resp_cert)[,2]
    resp_cert$pred <- resp_cert$var1.pred
    
    output_raster_cert <- terra::rasterize(resp_cert, raster_grid, field = "pred", fun = "mean")
    
    tm_shape(output_raster,
             bbox = st_bbox(output_raster)) + 
      tm_raster(col = "mean", col.scale = tm_scale_continuous(values = palette, limits = var_limits),
                col.legend = tm_legend(title = legend_title, show = legend_show, frame = FALSE, text.size = 1, title.size = 1, orientation = "landscape")) + 
      tm_layout(frame = FALSE,
                inner.margins = c(0,0,0,0)) + 
      tm_shape(output_raster_cert,
               bbox = st_bbox(output_raster)) + 
      tm_raster(col.scale = tm_scale_continuous(values = "grey"), col_alpha = "mean", col_alpha.scale = tm_scale_continuous(values = c(1,0)),
                col.legend = tm_legend(show = FALSE),
                col_alpha.legend = tm_legend(show = FALSE)) + 
      tm_layout(inner.margins = c(0,0,0,0),
                outer.margins = c(0,0,0,0)) + 
      tm_title(text = title_lbl, position = tm_pos_out("center","top", align.h = "left"), size = 1.2)
  }else{
    tm_shape(output_raster,
             bbox = st_bbox(output_raster)) + 
      tm_raster(col = "mean", col.scale = tm_scale_continuous(values = palette, limits = var_limits),
                col.legend = tm_legend(title = legend_title, show = legend_show, frame = FALSE, text.size = 1, title.size = 1, orientation = "landscape"))+
      tm_layout(frame = FALSE,
                inner.margins = c(0,0,0,0)) + 
      tm_title(text = title_lbl, position = tm_pos_out("center","top", align.h = "left"), size = 1.2)
  }
  
}
