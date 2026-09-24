source(here::here("R/bGEV/bGEVcode.R"))
source(here::here("data-analysis/R/predictive-utils.R"))
library(coda)

library(tidyverse)
library(sf)




STANDARD <- 35

## Some plotting object

california <- st_read(here::here("data-analysis/shp/ca_state/CA_State.shp"))
california <- st_transform(california, crs = 4326)
grid_clipped <- readRDS(here::here("data-analysis/results/grid_clipped.rds"))



#### PREFERENTIAL MODEL ####

mod1 <-
  readRDS(here::here("data-analysis/results/pm25_m1_thinned.rds"))

# generate using the mod1-posteriorpreds file
q_pred_mod1_2022 <- readRDS(here::here("data-analysis/results/q_pred_2022_mt1.rds"))
q_pred_mod1_2023 <- readRDS(here::here("data-analysis/results/q_pred_2023_mt1.rds"))
q_pred_mod1_2024 <- readRDS(here::here("data-analysis/results/q_pred_2024_mt1.rds"))


y_pred_mod1_2022 <- get_y_samples(mod1, q_pred_mod1_2022, 6)
y_pred_mod1_2023 <- get_y_samples(mod1, q_pred_mod1_2023, 7)
y_pred_mod1_2024 <- get_y_samples(mod1, q_pred_mod1_2024, 8)

saveRDS(y_pred_mod1_2024, here::here("data-analysis/results/y_pred_2024_mt1.rds"))

#### BASELINE MODEL ####

modbt1 <- readRDS(here::here("data-analysis/results/pm25_mb1_thinned.rds"))

q_pred_modbt1_2022 <- readRDS(here::here("data-analysis/results/q_pred_2022_mbt1.rds"))
q_pred_modbt1_2023 <- readRDS(here::here("data-analysis/results/q_pred_2023_mbt1.rds"))
q_pred_modbt1_2024 <- readRDS(here::here("data-analysis/results/q_pred_2024_mbt1.rds"))

y_pred_modbt1_2022 <- get_y_samples(modbt1, q_pred_modbt1_2022, 6)
y_pred_modbt1_2023 <- get_y_samples(modbt1, q_pred_modbt1_2023, 7)
y_pred_modbt1_2024 <- get_y_samples(modbt1, q_pred_modbt1_2024, 8)
saveRDS(y_pred_modbt1_2024, here::here("data-analysis/results/y_pred_2024_mbt1.rds"))

## Probability of Std Attainment

get_attainment_probs <- function(y1, y2, y3) {
  avg3 <- (exp(y1) + exp(y2) + exp(y3)) / 3
  probs <- apply(avg3, 2, function(x) mean(x > STANDARD))
  return(probs)
}

probs_m1_pred <- get_attainment_probs(y1 = y_pred_mod1_2022, y2 = y_pred_mod1_2023, y3 = y_pred_mod1_2024)
saveRDS(probs_m1_pred, here::here("data-analysis/results/probs_m1_pred.rds"))
probs_mb1_pred <- get_attainment_probs(y1 = y_pred_modbt1_2022, y2 = y_pred_modbt1_2023, y3 = y_pred_modbt1_2024)
saveRDS(probs_mb1_pred, here::here("data-analysis/results/probs_mb1_pred.rds"))
## MAPS

probs_m1_pred <- readRDS(here::here("data-analysis/results/probs_m1_pred.rds"))
probs_mb1_pred <- readRDS(here::here("data-analysis/results/probs_mb1_pred.rds"))

plot_sf <-
  st_as_sf(
    data.frame(
      "Baseline" = probs_mb1_pred,
      "Preferential" = probs_m1_pred,
      "Difference" = probs_m1_pred - probs_mb1_pred
    ),
    geometry = grid_clipped
  ) %>%
  pivot_longer(
    cols = c("Baseline", "Preferential", "Difference")
  ) %>%
  mutate(
    larger_70 = value > 0.8,
    larger_70 = value > 0.7,
    larger_50 = value > 0.5,
    larger_25 = value > 0.25
  )

# Make plot


probs_map <- ggplot() +
  geom_sf(
    data = california,
    fill = "gray99",
    alpha = 1,
    col = "black"
  ) +
  geom_sf(data = grid_clipped) +
  geom_sf(data = plot_sf %>% filter(name != "Difference"), aes(fill = value)) +
  facet_wrap(~name) +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    limits = c(0, 1),
    name = "Probability\nValue > 35 ug/m3"
  ) +
  theme_bw() +
  theme_bw(base_size = 14, base_family = "serif") + # use serif font
  theme(
    strip.text = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )


diff_map <- ggplot() +
  geom_sf(
    data = california,
    fill = "gray99",
    alpha = 1,
    col = "black"
  ) +
  geom_sf(data = grid_clipped) +
  geom_sf(data = plot_sf %>% filter(name == "Difference"), aes(fill = value)) +
  scale_fill_distiller(palette = "Spectral", name = "(Probs for Preferential - \nProbs for Baseline)") +
  theme_bw() +
  ggtitle("Difference in Probability of Non-Attainment")

pred_non_att_70 <- ggplot() +
  geom_sf(
    data = california,
    fill = "gray99",
    alpha = 1,
    col = "black"
  ) +
  geom_sf(data = grid_clipped) +
  geom_sf(data = plot_sf, aes(fill = larger_70)) +
  facet_wrap(~name) +
  scale_fill_grey(name = "Probability of \nNon-Attainment \n> 70%") +
  theme_bw()

pred_non_att_70 <- ggplot() +
  geom_sf(
    data = california,
    fill = "gray99",
    alpha = 1,
    col = "black"
  ) +
  geom_sf(data = grid_clipped) +
  geom_sf(data = plot_sf, aes(fill = larger_70)) +
  facet_wrap(~name) +
  scale_fill_grey(name = "Probability of \nNon-Attainment \n> 70%") +
  theme_bw()

pred_non_att_70 <- ggplot() +
  geom_sf(
    data = california,
    fill = "gray99",
    alpha = 1,
    col = "black"
  ) +
  geom_sf(data = grid_clipped) +
  geom_sf(data = plot_sf, aes(fill = larger_70)) +
  facet_wrap(~name) +
  scale_fill_grey(name = "Probability of \nNon-Attainment \n> 70%") +
  theme_bw()




#### For the Actual Stations
pm25 <- readRDS(here::here("data-analysis/pm2.5/clean/8years_sf.rds"))
dat10 <- readRDS(here::here("data-analysis/pm2.5/clean/10years_biggrid.rds"))

# grid
california <- st_read(here::here("data-analysis/shp/ca_state/CA_State.shp"))
california <- st_transform(california, crs = 4326)



##
grid <- st_make_grid(california, n = c(32, 32))
grid_clipped <- st_intersection(grid, california)
grid_centroids <- st_centroid(st_make_valid(grid_clipped))
pred_scaled_coords <- grid_centroids %>%
  st_coordinates() %>%
  data.frame() %>%
  mutate(
    x = (X - bbox_points[1]) * scale_x,
    y = (Y - bbox_points[2]) * scale_y
  ) %>%
  dplyr::select(x, y) %>%
  as.matrix()


grid_contains <- st_contains(st_make_valid(grid_clipped), pm25 %>% group_by(stations) %>% summarize(first(geometry)))

# Find intersections: returns a list, one element per point
intersection_list <- st_intersects(pm25 %>% group_by(stations) %>% summarize(first(geometry)), st_make_valid(grid_clipped))

# Extract the grid index for each point (assuming each point lies in at most one grid)
grid_index <- sapply(intersection_list, function(x) if (length(x) > 0) x[1] else NA)

# dat10 <- readRDS(here::here("data-analysis/pm2.5/clean/10years_biggrid.rds"))
tmp_df <- data.frame(station = dat10$stations, y = dat10$y, year = dat10$years) %>% filter(year %in% c(2022, 2023, 2024))

actual_std_att <- tmp_df %>%
  group_by(station) %>%
  summarize(mean_exp_y = mean(exp(y))) %>%
  mutate(passed = mean_exp_y > STANDARD)


shared_stations <- tmp_df %>%
  distinct(year, station) %>% # remove duplicate rows if they exist
  group_by(station) %>%
  summarize(n_years = n_distinct(year)) %>%
  filter(n_years == n_distinct(tmp_df$year)) %>%
  pull(station)

actual_std_att <- actual_std_att %>% filter(station %in% shared_stations)

names_keep <- paste0("q_", shared_stations)

q_est_modbt1_2022 <- readRDS(here::here("data-analysis/results/q_est_mbt1_2022.rds"))[, names_keep]
q_est_modbt1_2023 <- readRDS(here::here("data-analysis/results/q_est_mbt1_2023.rds"))[, names_keep]
q_est_modbt1_2024 <- readRDS(here::here("data-analysis/results/q_est_mbt1_2024.rds"))[, names_keep]

y_est_modbt1_2022 <- get_y_samples(modbt1, q_est_modbt1_2022, 6)
y_est_modbt1_2023 <- get_y_samples(modbt1, q_est_modbt1_2023, 7)
y_est_modbt1_2024 <- get_y_samples(modbt1, q_est_modbt1_2024, 8)


probs_mb1_est <- get_attainment_probs(y1 = y_est_modbt1_2022, y2 = y_est_modbt1_2023, y3 = y_est_modbt1_2024)
saveRDS(probs_mb1_est, here::here("data-analysis/results/probs_mb1_est.rds"))

q_est_mod1_2022 <- readRDS(here::here("data-analysis/results/q_est_2022_mt1.rds"))[, names_keep]
q_est_mod1_2023 <- readRDS(here::here("data-analysis/results/q_est_2023_mt1.rds"))[, names_keep]
q_est_mod1_2024 <- readRDS(here::here("data-analysis/results/q_est_2024_mt1.rds"))[, names_keep]

y_est_mod1_2022 <- get_y_samples(mod1, q_est_mod1_2022, 6)
y_est_mod1_2023 <- get_y_samples(mod1, q_est_mod1_2023, 7)
y_est_mod1_2024 <- get_y_samples(mod1, q_est_mod1_2024, 8)

probs_m1_est <- get_attainment_probs(y1 = y_est_mod1_2022, y2 = y_est_mod1_2023, y3 = y_est_mod1_2024)
# saveRDS(probs_m1_est, here::here("data-analysis/results/probs_m1_est.rds"))
probs_m1_est <- readRDS(here::here("data-analysis/results/probs_m1_est.rds"))
probs_mb1_est <- readRDS(here::here("data-analysis/results/probs_mb1_est.rds"))
# probs_mb1_est <- readRDS(here::here("data-analysis/results/probs_mb1_est.rds"))


## Add population
library(tidycensus)
ca_population <- get_acs(
  geography = "tract",
  state = "CA",
  variables = "B01003_001",
  year = 2023,
  geometry = T
)


ca_population <- st_transform(ca_population, crs = st_crs(pm25_sf))
stations_pop <- pm25_sf %>%
  filter(year %in% c(2022, 2023, 204)) %>%
  dplyr::select(stations) %>%
  st_join(ca_population %>% dplyr::select(estimate)) %>%
  distinct(stations, estimate) %>%
  rename(station = stations) %>%
  filter(station %in% shared_stations)


actual_std_att <- actual_std_att %>% left_join(stations_pop)

actual_std_att <- actual_std_att %>% mutate(
  base_25 = probs_mb1_est > 0.25,
  pref_25 = probs_m1_est > 0.25,
  base_70 = probs_mb1_est > 0.5,
  pref_70 = probs_m1_est > 0.5,
  base_70 = probs_mb1_est > 0.70,
  pref_70 = probs_m1_est > 0.70,
  estimate = estimate / 1000,
  base_80 = probs_mb1_est > 0.8,
  pref_80 = probs_m1_est > 0.8,
)


# tab1 <- actual_std_att %>% dplyr::select(passed, base_70) %>% group_by(base_70) %>% summarize(actual_exc = sum(passed), actual_not_exc = sum(passed==FALSE)) %>% arrange(desc(base_70)) %>% rename("predicted" = base_70) %>% mutate(model = "Baseline")  %>% dplyr::select(model, predicted, actual_exc, actual_not_exc)

tab1 <- actual_std_att %>%
  dplyr::select(passed, base_70, estimate) %>%
  group_by(base_70) %>%
  summarize(actual_exc = sum(passed), actual_not_exc = sum(passed == FALSE), pop_exc = sum(estimate[passed == T]), pop_not_exc = sum(estimate[passed == F])) %>%
  arrange(desc(base_70)) %>%
  rename("predicted" = base_70) %>%
  mutate(model = "Baseline") %>%
  dplyr::select(model, predicted, actual_exc, actual_not_exc, pop_exc, pop_not_exc)


# tab1_pop <- actual_std_att %>% dplyr::select(passed, base_70,estimate) %>% group_by(base_70,passed) %>% summarize(pop_exc = sum(estimate[passed==T]), pop_not_exc = sum(estimate[passed==F]))


tab2 <- actual_std_att %>%
  dplyr::select(passed, pref_70, estimate) %>%
  group_by(pref_70) %>%
  summarize(actual_exc = sum(passed), actual_not_exc = sum(passed == FALSE), pop_exc = sum(estimate[passed == T]), pop_not_exc = sum(estimate[passed == F])) %>%
  arrange(desc(pref_70)) %>%
  rename("predicted" = pref_70) %>%
  mutate(model = "Preferential") %>%
  dplyr::select(model, predicted, actual_exc, actual_not_exc, pop_exc, pop_not_exc)


## make table
tab <- rbind(tab1, tab2)
library(kableExtra)
tab %>%
  mutate(
    actual_exc = sprintf("%d (%.2f)", actual_exc, pop_exc),
    actual_not_exc = sprintf("%d (%.2f)", actual_not_exc, pop_not_exc)
  ) %>%
  dplyr::select(-pop_exc, -pop_not_exc) %>%
  rename(
    "Model" = model,
    "Prob Non-Att. > 70\\%" = predicted,
    "Std Exceeded" = actual_exc,
    "Std Not Exceeded" = actual_not_exc
  ) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    align = "lccc",
    escape = FALSE,
    caption = "Observed and Predicted NAAQS Std Attainment, Number of Stations (Census tract 2023 Population (in Thousands))",
    label = "tab::naaqstab1",
    col.names = c(
      "\\textbf{Model}",
      "\\textbf{Prob Non-Att.>0\\%}",
      "\\textbf{Std Exceeded}",
      "\\textbf{Std Not Exceeded}"
    )
  ) %>%
  kable_styling(latex_options = c("hold_position", "repeat_header", position = "bottom")) %>%
  collapse_rows(columns = 1, latex_hline = "none", valign = "middle")



true_pos_m1 <- c()
false_pos_m1 <- c()
true_pos_mb1 <- c()
false_pos_mb1 <- c()


probs_check <- seq(0.0, 1.0, 0.05)


pos_id <- which(actual_std_att$passed == TRUE)
i <- 0
for (k in probs_check) {
  i <- i + 1
  true_pos_m1[i] <- sum(((probs_m1_est > k) == actual_std_att$passed)[pos_id]) / length(actual_std_att$passed[pos_id])
  true_pos_mb1[i] <- sum(((probs_mb1_est > k) == actual_std_att$passed)[pos_id]) / length(actual_std_att$passed[pos_id])
  false_pos_m1[i] <- sum(((probs_m1_est > k) != actual_std_att$passed)[pos_id]) / length(actual_std_att$passed[pos_id])
  false_pos_mb1[i] <- sum(((probs_mb1_est > k) != actual_std_att$passed)[pos_id]) / length(actual_std_att$passed[pos_id])
}


plot(probs_check, true_pos_m1, "l")
lines(probs_check, true_pos_mb1)


plot_df <- rbind(
  data.frame("true_pos" = true_pos_m1, "false_pos" = false_pos_m1, model = "Preferential", "prob" = probs_check),
  data.frame("true_pos" = true_pos_mb1, "false_pos" = false_pos_mb1, model = "Baseline", "prob" = probs_check)
) %>% pivot_longer(cols = c("true_pos", "false_pos"))

plot_df$name <- factor(plot_df$name, levels = c("true_pos", "false_pos"))



ggplot(plot_df) +
  geom_line(aes(x = prob, y = value, color = model, group = model)) +
  geom_point(aes(x = prob, y = value, color = model)) +
  facet_wrap(~name, labeller = labeller(name = c(true_pos = "True Positive Rate", false_pos = "False Positive Rate"))) +
  theme_bw() +
  xlab("Probability Threshold") +
  ylab("") +
  scale_color_manual(name = "", values = c("Preferential" = "#E69F00", "Baseline" = "#0072B2")) +
  ggtitle("Attainment Analysis Classification") +
  theme_bw(base_size = 14, base_family = "serif") + # use serif font
  theme(
    strip.text = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )
