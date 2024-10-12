
# Archive Take Home Ex 01

# 2. Spatial Point Pattern Analysis

Before any meaningful analysis between any factors and accidents can be made, we need to know if the accidents events occur following a pattern or at random via spatial point patterns analysis. This analysis can also help us to understand if the events vary due to changes in underlying property (e.g. road types) or if they vary due to interaction effects between points. This will be established through first-order and second-order spatial point pattern analysis.

## 3.1 First-Order Spatial Point Pattern Analysis

Under this section, Kernel Density Estimation and Nearest Neighbour Distance will be explored as density-based and distance-based techniques respectively. This is bearing in mind their respective strengths.

### 3.1.1 Density-Based: Kernel Density Estimation

Additional preparation is required to perform the analysis. These include:\
- Preparing the accident data in ppp object for `spatstat`\
- Checking the accident data for duplicates\
- Preparing the owin object\

```{r}
# Set seed to ensure reproducibility of all subsequent analysis
set.seed(42)

# Convert accidents data to ppp format
ra_ppp <- as.ppp(ra_sf)

# Check accidents events for duplicates
any(duplicated(ra_ppp))
```

```{r}
# Create owin object
th_bound_owin <- as.owin(th_bound_sf)

# Combine events object and owin object
ra_th_bound_ppp <-  ra_ppp[th_bound_owin]

# Rescale ppp from meter to kilometer
ra_th_bound_ppp_km <- rescale.ppp(ra_th_bound_ppp, 1000, "km")
```

### 3.1.1.1 Adaptive Kernel Density Estimation

Based on the context of road accidents, adaptive kernel density estimation is more appropriate than fixed kernel density estimation as the accidents are more heterogeneous as observed from the EDA.

```{r}
# Perform adative bandwith KDE
kde_ra_adp <- adaptive.density(ra_th_bound_ppp_km, method="kernel")
plot(kde_ra_adp, 
     main = paste0("Adaptive KDE"),
     cex.axis = 0.8)
```

The intensities can be observed as small strips, which are notable roads (i.e. Kanchanaphisek Road) where sizable proportions of the accidents happen. By identifying these hotspots, further investigation and policies can be made to mitigate them.

Lets see if there are any differences when looking at individual provinces.

::: panel-tabset
# Bangkok

```{r}
# Create Bangkok ppp object
ra_bk_ppp <- ra_sf %>%
  filter(province_en == 'Bangkok') %>%
  as.ppp()

# Create Bangkok owin object
bk_owin <- th_bound_sf %>%
  filter(ADM1_EN == 'Bangkok') %>%
  as.owin()

# Combine events object and owin object
ra_bk_bound <-  ra_bk_ppp[bk_owin]

# Rescale ppp from meter to kilometer
ra_bk_bound_km <- rescale.ppp(ra_bk_bound, 1000, "km")

# Perform adaptive bandwidth KDE
kde_ra_adp_bk <- adaptive.density(ra_bk_bound_km, method="kernel")
plot(kde_ra_adp_bk, 
     main = paste0("Adaptive KDE: Bangkok"),
     cex.axis = 0.8)

```

# Nonthaburi

```{r}
# Create Nonthaburi ppp object
ra_ntb_ppp <- ra_sf %>%
  filter(province_en == 'Nonthaburi') %>%
  as.ppp()

# Create Nonthaburi owin object
ntb_owin <- th_bound_sf %>%
  filter(ADM1_EN == 'Nonthaburi') %>%
  as.owin()

# Combine events object and owin object
ra_ntb_bound <-  ra_ntb_ppp[ntb_owin]

# Rescale ppp from meter to kilometer
ra_ntb_bound_km <- rescale.ppp(ra_ntb_bound, 1000, "km")

# Perform adaptive bandwidth KDE
kde_ra_adp_ntb <- adaptive.density(ra_ntb_bound_km, method="kernel")
plot(kde_ra_adp_ntb, 
     main = paste0("Adaptive KDE: Nonthaburi"),
     cex.axis = 0.8)

```

# Pathum Thani

```{r}
# Create Pathum Thani ppp object
ra_pt_ppp <- ra_sf %>%
  filter(province_en == 'Pathum Thani') %>%
  as.ppp()

# Create Pathum Thani owin object
pt_owin <- th_bound_sf %>%
  filter(ADM1_EN == 'Pathum Thani') %>%
  as.owin()

# Combine events object and owin object
ra_pt_bound <-  ra_pt_ppp[pt_owin]

# Rescale ppp from meter to kilometer
ra_pt_bound_km <- rescale.ppp(ra_pt_bound, 1000, "km")

# Perform adaptive bandwidth KDE
kde_ra_adp_pt <- adaptive.density(ra_pt_bound_km, method="kernel")
plot(kde_ra_adp_pt, 
     main = paste0("Adaptive KDE: Pathum Thani"),
     cex.axis = 0.8)

```

# Samut Prakan

```{r}
# Create Samut Prakan ppp object
ra_sp_ppp <- ra_sf %>%
  filter(province_en == 'Samut Prakan') %>%
  as.ppp()

# Create Samut Prakan owin object
sp_owin <- th_bound_sf %>%
  filter(ADM1_EN == 'Samut Prakan') %>%
  as.owin()

# Combine events object and owin object
ra_sp_bound <-  ra_sp_ppp[sp_owin]

# Rescale ppp from meter to kilometer
ra_sp_bound_km <- rescale.ppp(ra_sp_bound, 1000, "km")

# Perform adaptive bandwidth KDE
kde_ra_adp_sp <- adaptive.density(ra_sp_bound_km, method="kernel")
plot(kde_ra_adp_sp, 
     main = paste0("Adaptive KDE: Samut Prakan"),
     cex.axis = 0.8)

```

# Samut Sakhon

```{r}
# Create Samut Sakhon ppp object
ra_ss_ppp <- ra_sf %>%
  filter(province_en == 'Samut Sakhon') %>%
  as.ppp()

# Create Samut Sakhon owin object
ss_owin <- th_bound_sf %>%
  filter(ADM1_EN == 'Samut Sakhon') %>%
  as.owin()

# Combine events object and owin object
ra_ss_bound <-  ra_ss_ppp[ss_owin]

# Rescale ppp from meter to kilometer
ra_ss_bound_km <- rescale.ppp(ra_ss_bound, 1000, "km")

# Perform adaptive bandwidth KDE
kde_ra_adp_ss <- adaptive.density(ra_ss_bound_km, method="kernel")
plot(kde_ra_adp_ss, 
     main = paste0("Adaptive KDE: Samut Sakhon"),
     cex.axis = 0.8)

```

# Nakhon Pathom

```{r}
# Create Nakhon Pathom ppp object
ra_np_ppp <- ra_sf %>%
  filter(province_en == 'Nakhon Pathom') %>%
  as.ppp()

# Create Nakhon Pathom owin object
np_owin <- th_bound_sf %>%
  filter(ADM1_EN == 'Nakhon Pathom') %>%
  as.owin()

# Combine events object and owin object
ra_np_bound <-  ra_np_ppp[np_owin]

# Rescale ppp from meter to kilometer
ra_np_bound_km <- rescale.ppp(ra_np_bound, 1000, "km")

# Perform adaptive bandwidth KDE
kde_ra_adp_np <- adaptive.density(ra_np_bound_km, method="kernel")
plot(kde_ra_adp_np, 
     main = paste0("Adaptive KDE: Nakhon Pathom"),
     cex.axis = 0.8)

```
:::
  
  Conclusion:
  
  ### 3.1.1.2 Fixed Kernel Density Estimation
  
  While fixed KDE is less applicable, it will be interesting to compare the output with adaptive KDE.

`bw.ppl()` will be used as the automatic kernel method to determine the sigma as the pattern seems to consist mostly of tight clusters around major roads.\
Gaussian will be the smoothing kernel of choice, as it is widely used and there are no negative events in the current context.

```{r}
# Get fixed bandwidth 
fixed_bw <- bw.ppl(ra_th_bound_ppp_km)

# Perform fixed kernel density estimation
plot(density(ra_th_bound_ppp_km, sigma=fixed_bw, edge=TRUE, kernel="gaussian"), 
     main = paste0("Fixed KDE: Bandwidth=", round(fixed_bw, 4)),
     cex.axis = 0.8)
```

Interestingly,

### 3.1.2 Distance-Based: Nearest Neighbour Index

`clarkevans.test()` of spatstat will be used to perform the Clark-Evans test of aggregation.\
From the previous observations of the road accidents, it is likely that there are clustering involved, therefore, the alternative for hypothesis test is "clustered".

H0: Distribution of road accidents are randomly distributed\
H1: Distribution of road accidents are *not* randomly distributed (corresponding to a clustered point pattern)

```{r}
# Perform Clark-Evans on study area 
clarkevans.test(ra_th_bound_ppp,
                correction="none",
                clipregion="th_bound_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```

Conclusion: As R is significantly smaller than 1 (0.19), this indicates strong clustering. p-value being \< 0.05 allows us to reject the null hypothesis and accept the alternative hypothesis that the point pattern is clustered and not random.

Now lets do the same for each province

::: panel-tabset
# Bangkok

```{r}
# Perform Clark-Evans on Bangkok
clarkevans.test(ra_bk_bound,
                correction="none",
                clipregion="bk_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```

# Nonthaburi

```{r}
# Perform Clark-Evans on Nonthaburi
clarkevans.test(ra_ntb_bound,
                correction="none",
                clipregion="ntb_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```

# Pathum Thani

```{r}
# Perform Clark-Evans on Pathum Thani
clarkevans.test(ra_pt_bound,
                correction="none",
                clipregion="pt_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```

# Samut Prakan

```{r}
# Perform Clark-Evans on Samut Prakan
clarkevans.test(ra_sp_bound,
                correction="none",
                clipregion="sp_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```

# Samut Sakhon

```{r}
# Perform Clark-Evans on Samut Sakhon
clarkevans.test(ra_ss_bound,
                correction="none",
                clipregion="ss_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```

# Nakhon Pathom

```{r}
# Perform Clark-Evans on Nakhon Pathom
clarkevans.test(ra_np_bound,
                correction="none",
                clipregion="np_owin",
                alternative=c("clustered"),
                nsim=99 # 100 simulation runs
)
```
:::
  
  Conclusion: Although the R values for the provinces are generally larger than that for the entire study area, they are still much smaller than 1, indicating strong clustering. p-values for all provinces being \< 0.05 allows us to reject the null hypothesis and accept the alternative hypothesis that the point pattern is clustered and not random within each province.

## 3.2 Second-Order Spatial Point Pattern Analysis

G, F and L functions will be explored as they can each give useful insights:\
- G function (*measures the distribution of the distances from an arbitrary event to its nearest event*): Detect local clustering of accidents to identify broader underlying hazards.
- L function (: Understand clustering of accidents across multiple spatial scales
              
              
              
              ### 3.2.1 G Function
              
              ```{r}
              # Compute G function
              G_th = Gest(ra_th_bound_ppp, correction = "best")
              plot(G_th, xlim=c(0,500))
              ```
              
              ```{r}
              # Complete spatial randomness test with G function
              G_th.csr <- envelope(ra_th_bound_ppp, Gest, nsim = 99) # 100 simulations 
              plot(G_th.csr)
              ```
              
              
              ### 3.2.2 F Function
              
              ```{r}
              # Compute F function
              F_th = Fest(ra_th_bound_ppp, correction = "best")
              plot(F_th)
              ```
              
              ```{r}
              # Complete spatial randomness test with G function
              F_th.csr <- envelope(ra_th_bound_ppp, Fest, nsim = 99) # 100 simulations 
              plot(F_th.csr)
              ```
              
              
              ## 3.3 Conclusion
              
              
              
              
              
              
              
              
              
              
              
              ```{r}
              #| eval: false # eval: false
              future::plan(future::multisession(workers=8))
              
              # Create lixels
              th_lixels_bk <- lixelize_lines(bk_osm_sf, 200)
              
              # Create samples
              th_samples_bk <- lines_center(th_lixels_bk) 
              
              ra_bk_sf <- ra_sf %>%
                filter(province_en == "Bangkok")
              
              th_densities_bk <- nkde.mc(bk_osm_sf, 
                                         events = ra_bk_sf,
                                         w = rep(1, nrow(ra_bk_sf)),
                                         samples = th_samples_bk,
                                         kernel_name = "quartic",
                                         bw = 250, 
                                         div= "bw", 
                                         method = "simple", 
                                         digits = 1, 
                                         tol = 0.1,
                                         grid_shape = c(5,5), 
                                         max_depth = 8,
                                         agg = 5, 
                                         sparse = TRUE, # Slower but require less memory
                                         verbose = FALSE)
              ```
              
              ```{r}
              #| eval: false # eval: false
              th_samples_bk$density <- th_densities_bk
              th_lixels_bk$density <- th_densities_bk
              
              th_samples_bk$density <- th_samples_bk$density*1000
              th_lixels_bk$density <- th_lixels_bk$density*1000
              ```
              
              ```{r, fig.width=12, fig.height=8}
              #| eval: false # eval: false
              my_palette <- c("blue", "green")
              tmap_mode('plot')
              tm_shape(th_lixels_bk)+
                tm_lines(col="density", palette = my_palette, style = "cont") +
                tm_shape(ra_bk_sf)+
                tm_dots(size = 0.001, col = 'yellow')
              ```
              
              ```{r}
              library(leaflet)
              library(crosstalk)
              library(sf)
              library(dplyr)
              library(RColorBrewer)
              
              # Ensure 'highway' column is a factor
              th_osm_sf$highway <- as.factor(th_osm_sf$highway)
              
              th_osm_sf_new <- th_osm_sf %>%
                st_transform(crs = 4326)
              
              # Convert the road sf dataframe to SharedData for filtering
              shared_road_sf <- SharedData$new(th_osm_sf_new)
              
              # Define a color palette for the road types (highway column)
              pal_highway <- colorFactor(brewer.pal(n = 6, name = "Set1"), domain = th_osm_sf$highway)
              
              # Create the interactive map with filters using crosstalk and leaflet
              bscols(
                widths = c(3, 9),
                # Add the checkbox filter for "highway"
                list(
                  filter_checkbox("highway", "Road Type", shared_road_sf, ~highway, columns = 1)
                ),
                
                # Add leaflet map with the filtered data
                leaflet(shared_road_sf) %>% 
                  addProviderTiles(providers$CartoDB.Positron) %>%
                  # Plot the roads, colored by the "highway" column
                  addPolylines(color = ~pal_highway(highway), weight = 2, opacity = 0.7) 
              )
              ```
              
              
              library(spgwr)
              
              # Example: GWR with weather condition as a predictor and spatial coordinates
              gwr_model <- gwr(severe_injury ~ weather_condition + vehicle_type, 
                               data = accidents_sf, coords = cbind(accidents_sf$X, accidents_sf$Y), 
                               bandwidth = 50)
              
              summary(gwr_model)
              
              
              ```{r}
              # Assuming 'coords' contains the coordinates of your data points (province centroids)
              max_k <- 15  # Set a maximum value of k to test
              distances <- sapply(1:max_k, function(k) {
                mean(spatstat.geom::nndist(coords_cent, k = k))  # Calculate the average distance to the k-th neighbor
              })
              
              # Plot the elbow curve
              elbow_data <- data.frame(k = 1:max_k, distance = distances)
              ggplot2::ggplot(elbow_data, aes(x = k, y = distance)) +
                geom_line() +
                geom_point() +
                labs(title = "Elbow Method for Determining k",
                     x = "Number of Neighbors (k)",
                     y = "Average Distance to k-th Neighbor")
              ```

              
              
              ```{r}
              # Step 1: Find the centroid of Phuket (or use the entire geometry)
              phuket_geom <- th_bound_tour_all %>%
                filter(ADM1_EN == "Phuket") %>%
                st_geometry()
              
              # Step 2: Find the centroids of all other regions
              other_geoms <- th_bound_tour_all %>%
                filter(ADM1_EN != "Phuket") %>%
                st_geometry()
              
              ```
              
              ```{r}
              # Step 3: Calculate the nearest neighbor using distance
              nearest_region_index <- st_nearest_feature(phuket_geom, other_geoms)
              
              phuket_index <- which(th_bound_tour_all$ADM1_EN == "Phuket")
              nb[[phuket_index]] <- c(nearest_region_index)
              nb[[nearest_region_index]] <- c(nb[[nearest_region_index]], phuket_index)
              ```
              
              
              # Function to 
              local_moran_space <- function(data, col, value, kpi) {
                # Set Seed
                set.seed(42)
                
                # Filter sf dataframe
                filtered_data <- data %>%
                  filter(!!sym(col) == value)
                
                lisa_ogp <- th_tour_knn5_all %>%
                  mutate(local_moran = local_moran(
                    tourist_stay_all, nb, wt, nsim = 999),
                    .before = 1) %>%
                  unnest(local_moran)
                
                
                # Append neighbour list and calculate spatial weight matrix
                filtered_data <- filtered_data %>%
                  mutate(nb = nb_knn_5,
                         wt = st_weights(nb, style = "W"),
                         wt_inv = st_inverse_distance(include_self(nb), geometry = coords_cent, 
                                                      scale = 100, alpha = 1),
                         local_moran = local_moran(!!sym(kpi), nb, wt, nsim = 999),
                         .before = 1) %>%
                  unnest(local_moran)
                
                # map1 <- tm_shape(filtered_data) +
                #   tm_fill("ii") + 
                #   tm_borders(alpha = 0.5) +
                #   tm_view(set.zoom.limits = c(6,8)) +
                #   tm_layout(main.title = paste0("local Moran's I ",kpi," (",value, ")"),
                #             main.title.size = 0.8)
                # 
                # map2 <- tm_shape(filtered_data) +
                #   tm_fill("p_ii",
                #         breaks = c(0, 0.001, 0.01, 0.05, 1),
                #             labels = c("0.001", "0.01", "0.05", "Not sig")) + 
                #   tm_borders(alpha = 0.5) +
                #   tm_layout(main.title = paste0("p-value of local Moran's I ",kpi," (",value, ")"),
                #             main.title.size = 0.8)
                
                filter_sig <- filtered_data  %>%
                  filter(p_ii_sim < 0.05)
                
                map3 <- tm_shape(filtered_data) +
                  tm_polygons() +
                  tm_borders(alpha = 0.5) +
                  tm_shape(filter_sig) +
                  tm_fill("mean") + 
                  tm_borders(alpha = 0.4) + 
                  tm_layout(main.title = paste0("LISA Map of ",kpi," (",value, ")"),
                            main.title.size = 0.4,
                            legend.text.size = 0.3,
                            legend.title.size = 0.4)
                # tm_layout(main.title = paste0("LISA Map of ",kpi," (",value, ")"),
                #       main.title.size = 0.8)
                
                # all_map <- tmap_arrange(map1, map2, map3, ncol = 3)
                # return(all_map)
                # return(list(map1 = map1, map2 = map2, map3 = map3))
                return(map3)