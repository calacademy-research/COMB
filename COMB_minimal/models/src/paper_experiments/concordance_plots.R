library(purrr)
library(dplyr)
library(ltm)

jags_results_dir <- "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_single_covar_all_sp/"
date_patterns <- c("2025-08-25.RData","2025-08-26.RData","2025-08-27.RData", "2025-08-28.RData")
regex_pattern <- paste0("(", paste(date_patterns, collapse = "|"), ")$")
my_jags_files <- list.files(path = jags_results_dir, pattern = regex_pattern, full.names = FALSE)

quantile(jagsResult$sims.list$mu[,2] - jagsResult$sims.list$mu[,1], probs = c(0.025,0.975)) 

get_components_from_filename <- function(filename, 
                                         sp_name_index = 3, 
                                         model_index = 6, 
                                         th_index = 2, 
                                         year_index = 4){
  components <- str_split(filename, pattern = "_")
  
  return(list(sp_name = components[[1]][sp_name_index], 
              model_name = paste("model", components[[1]][model_index], sep="_"), 
              threshold_val = components[[1]][th_index],
              year = components[[1]][year_index]
              ))
}

jags_results_list <- list()

finished_species <- names(jags_results_list)[-1]
for (filename in my_jags_files){
  metadata = get_components_from_filename(filename)
  species_name <- metadata$sp_name
  if (species_name %in% finished_species) {
    next
  }
  year_val <- metadata$year
  model_name <- metadata$model_name
  threshold_val <- metadata$threshold_val
  load(file.path(jags_results_dir,filename))
  jagsResult$sims.list <- NULL
  jagsResult$samples <- list()
  jags_results_list[[species_name]][[threshold_val]][[model_name]] <- jagsResult
  rm(jagsResult)
}

kos <- c("model_A", "model_AS", "model_H", "model_HA",
         "model_HAS", "model_HS", "model_S")
species_set <- names(jags_results_list)
thresholds <- names(jags_results_list[[1]])


for (sp in species_set){
  for (th in thresholds){
    for (ko in kos){ 
      jags_results_list[[sp]][[th]][[ko]]$samples <- list()
}}}


names(jags_results_list)
names(jags_results_list[["BBWO"]])
names(jags_results_list[["BBWO"]][["-2"]])
names(jags_results_list[["BBWO"]][["-2"]][["model_HAS"]])
rownames(jags_results_list[["BBWO"]][["-2"]][["model_HAS"]]$summary)


jags_beta1_df <- data.frame()

for (sp in species_set){
  for (th in thresholds){
    for (ko in kos){
      if ("beta1" %in% names(jags_results_list[[sp]][[th]][[ko]]$mean)){
      meta_df <- data.frame(species = sp, threshold = as.numeric(th), model = ko)
      beta1_results_df <- data.frame(t(jags_results_list[[sp]][[th]][[ko]]$summary["beta1",]))
      beta1_results_df$direction <- ifelse(beta1_results_df$mean > 0, 1, -1)
      beta1_results_df$overlap0direction <- beta1_results_df$direction * (1 - beta1_results_df$overlap0)
      jags_beta1_df  <- rbind.data.frame(jags_beta1_df, cbind.data.frame(meta_df, beta1_results_df))
    }}
  }
}

jags_det_param_df <- data.frame()
param_set <- c("mean_psi","p11","p_aru11","p_aru01","mu[1]","mu[2]", "beta1")
for (sp in species_set){
  for (th in thresholds){
    for (ko in kos){
      meta_df <- data.frame(species = sp, threshold = as.numeric(th), model = ko)
      for (param_name in param_set){
        if (param_name %in% rownames(jags_results_list[[sp]][[th]][[ko]]$summary)){
          det_param_df  <- data.frame(t(jags_results_list[[sp]][[th]][[ko]]$summary[param_name,]))
          det_param_df$param <- param_name
          jags_det_param_df   <- rbind.data.frame(jags_det_param_df, cbind.data.frame(meta_df, det_param_df))
        }
      }
    }
  }
}

jags_det_param_df %>% filter(param == "p_aru11")
jags_det_param_df %>% filter(species %in% target_species_set, Rhat > 1.2, threshold == -2)
## separate plots for threshold
## plot by groups of mean_psi
jags_det_param_df %>% filter(species %in% target_species_set, 
                             param == "beta1", threshold == -2,
                             Rhat > 1.2)


jags_det_param_df %>% 
  filter(model == "model_HA", param %in% c("p_aru11","p_aru01","p11")) %>%
  dplyr::select(species, model, threshold, mean, param) %>% 
  pivot_wider(id_cols = c("species", "model", "threshold"),
    names_from = "param",
              values_from = "mean") %>% 
  ggplot(aes(x=p11,y=p_aru11, label = species)) +
 # geom_point(aes(x=p11,y=p_aru11)) +
  geom_text(aes(label = species)) +
  coord_fixed(ratio = 1, xlim = c(0,0.8), ylim = c(0,0.8)) +
  facet_wrap(~threshold)

jags_det_param_df %>% filter(model == "model_HAS", param %in% c("mean_psi")) %>% group_by(threshold) %>% summarize(
  mean_mean_psi = mean(mean),
  quant25_mean_psi = quantile(mean, probs = 0.25),
  quant50_mean_psi = quantile(mean, probs = 0.50),
  quant75_mean_psi = quantile(mean, probs = 0.75)
)

jags_det_param_df %>% filter(param %in% c("mean_psi")) %>%
  ggplot(aes(y=mean, x = model, fill = as.factor(threshold))) +
  geom_boxplot()

jags_det_param_df %>% filter(param %in% c("mean_psi"), model == "model_HAS") %>%
  ggplot(aes(x=mean, fill = as.factor(threshold))) +
  geom_density()


jags_det_param_df %>% 
  filter(model == "model_HA", param %in% c("p_aru11","p_aru01","p11", "mean_psi")) %>%
  dplyr::select(species, model, threshold, mean, param) %>% 
  pivot_wider(id_cols = c("species", "model", "threshold"),
              names_from = "param",
              values_from = "mean") %>% 
  mutate(mean_psi_bucket = cut(mean_psi, 
                               breaks = c(0, 0.25, 0.75, 1.0), 
                               labels = c("Low: 0-0.25", "Medium: 0.25-0.75", "High: 0.75-1"),
                               right = FALSE) ) %>% 
  ggplot(aes(x=p11,y=p_aru11, label = species, color = mean_psi_bucket)) +
  # geom_point(aes(x=p11,y=p_aru11)) +
  geom_text(aes(label = species, color = mean_psi_bucket), size = 3, angle = 45, position = "nudge") +
  coord_fixed(ratio = 1, xlim = c(0,0.8), ylim = c(0,0.8)) +
  facet_wrap(~threshold)

target_species_set <- c( "BBWO", "CONI", "HETH", "HEWA", "MOUQ", "OSFL", "RBSA", "SOGR", "RBNU", "WETA", "LAZB", "BTYW")
REMOVE <- c("DUFL", "AMRO", "WEWP")

det_scatter_plot_2 <- jags_det_param_df %>% 
  filter(model == "model_HA", param %in% c("p_aru11","p_aru01","p11", "mean_psi"), 
         threshold == -2) %>%
  dplyr::select(species, model, threshold, mean, param) %>% 
  pivot_wider(id_cols = c("species", "model", "threshold"),
              names_from = "param",
              values_from = "mean") %>% 
  filter(mean_psi >= 0.05) %>%
  mutate(mean_psi_bucket = cut(mean_psi, 
                               breaks = c(0, 0.25, 0.75, 1.0), 
                               labels = c("Low: 0.05-0.25", "Medium: 0.25-0.75", "High: 0.75-1"),
                               right = FALSE) ) %>% 
  ggplot(aes(x=p11,y=p_aru11, label = species, color = mean_psi_bucket)) +
  labs(title = "Detection probability by method",
       subtitle = "Model = HAS, Threshold = -2",
       x = "Probability of Human Detection",
       y = "Probability of ARU detection") + 
  # geom_point(aes(x=p11,y=p_aru11)) +
  geom_text(aes(label = species, color = mean_psi_bucket), size = 3, angle = 45, position = "nudge")

ggsave("detect_probs_thresh-2.png", 
       det_scatter_plot_2,
       width = 6, height = 6, units = "in", dpi = 300)


det_scatter_plot_2a <- jags_det_param_df %>% 
  filter(model == "model_HA", param %in% c("p_aru11","p_aru01","p11", "mean_psi"), 
         threshold == -2, species %in% target_species_set) %>%
  dplyr::select(species, model, threshold, mean, param) %>% 
  pivot_wider(id_cols = c("species", "model", "threshold"),
              names_from = "param",
              values_from = "mean") %>% 
 # filter(mean_psi >= 0.05) %>%
  mutate(mean_psi_bucket = cut(mean_psi, 
                               breaks = c(0, 0.25, 0.75, 1.0), 
                               labels = c("Low: 0.0-0.25", "Medium: 0.25-0.75", "High: 0.75-1"),
                               right = FALSE) ) %>% 
  ggplot(aes(x=p11,y=p_aru11, label = species, color = mean_psi_bucket)) +
  labs(title = "Detection probability by method",
       subtitle = "Model = HAS, Threshold = -2, unfiltered",
       x = "Probability of Human Detection",
       y = "Probability of ARU detection") + 
  # geom_point(aes(x=p11,y=p_aru11)) +
  geom_text(aes(label = species, color = mean_psi_bucket), size = 3, angle = 45, position = "nudge")

ggsave("detect_probs_thresh-2a.png", 
       det_scatter_plot_2a,
       width = 6, height = 6, units = "in", dpi = 300)

det_scatter_plot_0 <- 
  jags_det_param_df %>% 
  filter(model == "model_HA", param %in% c("p_aru11","p_aru01","p11", "mean_psi"), threshold == 0) %>%
  dplyr::select(species, model, threshold, mean, param) %>% 
  pivot_wider(id_cols = c("species", "model", "threshold"),
              names_from = "param",
              values_from = "mean") %>% 
  filter(mean_psi >= 0.05) %>%
  mutate(mean_psi_bucket = cut(mean_psi, 
                               breaks = c(0, 0.25, 0.75, 1.0), 
                               labels = c("Low: 0.05-0.25", "Medium: 0.25-0.75", "High: 0.75-1"),
                               right = FALSE) ) %>% 
  ggplot(aes(x=p11,y=p_aru11, label = species, color = mean_psi_bucket)) +
  labs(title = "Detection probability by method",
       subtitle = "Model = HAS, Threshold = 0",
       x = "Probability of Human Detection",
       y = "Probability of ARU detection") + 
  # geom_point(aes(x=p11,y=p_aru11)) +
  geom_text(aes(label = species, color = mean_psi_bucket), size = 3, angle = 45, position = "nudge")
ggsave("detect_probs_thresh_0.png", 
       det_scatter_plot_0,
       width = 6, height = 6, units = "in", dpi = 300)






#### Count how frequently a model is ranked in each position within a species based on sd of beta1
count_sd_rank <- jags_beta1_df %>% 
  group_by(species) %>% 
  mutate(lowest_sd_rank = dense_rank(sd)) %>% 
  ungroup() %>%
  group_by(lowest_sd_rank, model) %>%
  tally() %>%
  pivot_wider(names_from = lowest_sd_rank, values_from = n) 

mean_rank_table <-  jags_beta1_df %>% 
  group_by(species) %>% 
  mutate(lowest_sd_rank = dense_rank(sd)) %>% 
  ungroup() %>%
  group_by(model) %>% 
  summarize(mean_rank = mean(lowest_sd_rank)) %>%
  arrange(mean_rank)

count_sd_rank %>% right_join(mean_rank_table, by="model") %>% arrange(mean_rank)

jags_beta1_df %>% 
  group_by(species) %>% 
  mutate(lowest_sd_rank = dense_rank(sd)) %>% ungroup() %>%
  group_by(lowest_sd_rank, model) %>%
  tally() %>%
  pivot_wider(names_from = lowest_sd_rank, values_from = n) 

mean_beta1_wide_df <- jags_beta1_df[,c("species","model","mean")] %>% 
  pivot_wider(
  id_cols = species,
  names_from=model,
  values_from=mean) %>%
  filter(!is.na(model_S))
mean_beta1_wide_mat <- base::as.matrix(mean_beta1_wide_df %>% column_to_rownames("species"))
overlap0_beta1_wide_df <- jags_beta1_df[,c("species","model","overlap0")] %>% pivot_wider(
  id_cols = species,
  names_from=model,
  values_from=overlap0) %>%
  filter(!is.na(model_S))
overlap0_beta1_wide_mat <- as.matrix(overlap0_beta1_wide_df %>% column_to_rownames("species"))

icc(mean_beta1_wide_mat, model = "oneway", type = "consistency")
icc(overlap0_beta1_wide_mat , model = "oneway", type = "agreement")

overlap0_beta1_wide_df %>%  mutate(pct_overlap0 = rowMeans(pick(where(is.numeric), -species)),
                                   cnt_overlap0 = rowSums(pick(where(is.numeric), 
                                                               -species, 
                                                               -pct_overlap0)))


library(ggpubr)
my_colors <- c("H" = "#59A14F",
               "HA" = "#76B7B2",
               "HS" = "#0072B2",
               "HAS" = "#CC79A7",
               "A" = "#E15759",
               "AS"="#FF9DA7",
               "S"="#E69F00")


jags_beta1_df <- jags_beta1_df %>% mutate(model_labs =str_split_i(model, "_", i = -1))
levels(jags_beta1_df$model_labs) <- c("H", "HA", "HS","HAS","A","AS","S")

beta1_plot_by_species <- 
jags_beta1_df %>% 
  # filter(threshold == 0) %>% 
  mutate(model_labs =factor(model_labs, levels=c("H", "HA", "HS","HAS","A","AS","S"))) %>%
  group_by(species) %>% mutate(mean_beta = mean(mean)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, mean_beta, .desc = TRUE)) %>%
  ggplot() +
  geom_point(aes(y=model_labs, x=`mean`, group = threshold, color=model_labs)) +
  geom_errorbar(aes(y=model_labs, xmin=X2.5., xmax=X97.5., color=model_labs,
                    linetype=as.factor(overlap0)), width=0.6) +
  scale_linetype_manual(values = c("0" = "dashed", "1" = "solid")) +
  facet_wrap(~species, strip.position = "right",ncol=1) +
 # geom_line(aes(x=mean, y=mean, color=modelLabs)) +
 # facet_wrap(~species) +
  theme_bw() +
  geom_vline(aes(xintercept = 0)) +
  labs(y="Model", 
       x = expression(
         paste(beta[1], ": relationship of fire severity with occupancy")),
       color = "Model",
       linetype = "HDI overlaps 0") +
  theme(legend.position = "right",
        legend.title.position = "top",
        # legend.byrow = TRUE,
        axis.text.y = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(size = 6),
        strip.text = element_text(face = "bold")) +
  guides(color = guide_legend(ncol=1)) +
  scale_color_manual(values = my_colors) 
ggsave("beta1_by_species.png", beta1_plot_by_species,
       width = 5, height = 10, units = "in", dpi = 300)

beta1_plot_by_all_species <- 
  jags_beta1_df %>% filter(threshold == -2)
  mutate(model_labs =factor(model_labs, levels=c("H", "HA", "HS","HAS","A","AS","S"))) %>%
  group_by(species) %>% mutate(mean_beta = mean(mean)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, mean_beta, .desc = TRUE)) %>%
  ggplot() +
  geom_point(aes(y=model_labs, x=`mean`, color=model_labs)) +
  geom_errorbar(aes(y=model_labs, xmin=X2.5., xmax=X97.5., color=model_labs,
                    linetype=as.factor(overlap0)), width=0.6) +
  scale_linetype_manual(values = c("0" = "dashed", "1" = "solid")) +
  facet_wrap(~species, strip.position = "right",ncol=3, dir = "v") +
  # geom_line(aes(x=mean, y=mean, color=modelLabs)) +
  # facet_wrap(~species) +
  theme_bw() +
  geom_vline(aes(xintercept = 0)) +
  labs(y="Model", 
       x = expression(
         paste(beta[1], ": relationship of fire severity with occupancy")),
       color = "Model",
       linetype = "HDI overlaps 0") +
  theme(legend.position = "right",
        legend.title.position = "top",
        # legend.byrow = TRUE,
        axis.text.y = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(size = 6),
        strip.text = element_text(face = "bold")) +
  guides(color = guide_legend(ncol=1)) +
  scale_color_manual(values = my_colors) 
ggsave("beta1_by_all_species.png", beta1_plot_by_all_species,
       width = 5, height = 10, units = "in", dpi = 300)


beta1_plot_by_target_species <- 
  jags_beta1_df %>% filter(threshold == -2, species %in% target_species_set) %>%
  mutate(model_labs =factor(model_labs, levels=c("H", "HA", "HS","HAS","A","AS","S"))) %>%
  group_by(species) %>% mutate(mean_beta = mean(mean)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, mean_beta, .desc = TRUE)) %>%
  ggplot() +
  geom_point(aes(y=fct_rev(model_labs), x=`mean`, color=model_labs)) +
  geom_errorbar(aes(y=fct_rev(model_labs), xmin=X2.5., xmax=X97.5., color=model_labs,
                    linetype=as.factor(overlap0)), width=0.6) +
  scale_linetype_manual(values = c("0" = "1242", "1" = "solid")) +
  facet_wrap(~species, strip.position = "right",ncol=3, dir = "v") +
  # geom_line(aes(x=mean, y=mean, color=modelLabs)) +
  # facet_wrap(~species) +
  theme_bw() +
  geom_vline(aes(xintercept = 0)) +
  labs(y="Model", 
       x = expression(
         paste("Posterior Mean and 95% HDI of ", beta[1], ": relationship of fire severity with occupancy")),
       color = "Model",
       linetype = "HDI overlaps 0") +
  theme(legend.position = "right",
        legend.title.position = "top",
        # legend.byrow = TRUE,
        axis.text.y = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(size = 6),
        strip.text = element_text(face = "bold")) +
  guides(color = guide_legend(ncol=1)) +
  scale_color_manual(values = my_colors) 
ggsave("beta1_by_target_species.png", beta1_plot_by_target_species,
       width = 8, height = 8, units = "in", dpi = 300)

my_colors <- c("H" = "#59A14F",
               "HA" = "#76B7B2",
               "HS" = "#0072B2",
               "HAS" = "#CC79A7",
               "A" = "#E15759",
               "AS"="#FF9DA7",
               "S"="#E69F00")

beta1_plot_by_model_column <- 
jags_beta1_df %>% mutate(model_labs =str_split_i(model, "_", i = -1)) %>%
  mutate(species_ordered = 
           reorder_within(species, 
                          as.numeric(mean), 
                          model_labs, 
                          fun = "identity")) %>%
  ggplot() +
  scale_y_reordered() +
  
  geom_point(aes(y=species_ordered, x=`mean`, color=species)) +
  geom_errorbar(aes(y=species_ordered, xmin=X2.5., xmax=X97.5., color=species,
                    linetype=as.factor(overlap0)), linewidth=.6
                    ) +
  scale_linetype_manual(values = c("0" = "dashed", "1" = "solid"))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title.position = "left",
        legend.direction = "horizontal",
        # legend.byrow = TRUE,
        axis.text.y = element_text(angle = 90,hjust=0.5,size=8),
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(size = 6)) +
  facet_wrap(~model_labs, strip.position = "top",nrow=1,scales = "free_y") +
  # geom_line(aes(x=mean, y=mean, color=modelLabs)) +
  # facet_wrap(~species) +
  geom_vline(aes(xintercept = 0)) +
  labs(y="Species code", 
       x = expression(
         paste(beta[1], ": relationship of fire severity with occupancy")),
       color="Species code",
       linetype = "HDI overlaps 0") +
  guides(color = guide_legend(nrow=1)) +
  theme(strip.text = element_text(face = "bold"))
ggsave("beta1_by_model_column.png", beta1_plot_by_model_column,
       width = 12, height = 6, units = "in", dpi = 300)

beta1_plot_by_model_row <- 
jags_beta1_df %>% mutate(model_labs =str_split_i(model, "_", i = -1)) %>%
  mutate(species_ordered = 
           reorder_within(species, 
                          as.numeric(mean), 
                          model_labs, 
                          fun = "identity")) %>%
  ggplot() +
  scale_y_reordered() +
  geom_point(aes(y=species_ordered, x=`mean`, color=species)) +
  geom_errorbar(aes(y=species_ordered, 
                    xmin=X2.5., 
                    xmax=X97.5., 
                    color=species,
                    linetype=as.factor(overlap0)), linewidth=.6) +
  scale_linetype_manual(values = c("0" = "dashed", "1" = "solid"))+
  facet_wrap(~model_labs, strip.position = "right",ncol=1,scales = "free_y") +
  # geom_line(aes(x=mean, y=mean, color=modelLabs)) +
  # facet_wrap(~species) +
  theme_bw() +
  geom_vline(aes(xintercept = 0)) +
  labs(y="Species", 
       x = expression(
         paste(beta[1], ": relationship of fire severity with occupancy")),
       color = "Species",
       linetype = "HDI overlaps 0") +
  theme(legend.position = "right",
        axis.text.y = element_text(angle = 0,size=8),
        strip.text = element_text(face = "bold") )+
  guides(color = guide_legend(ncol=1)) 
ggsave("beta1_by_model_row.png", beta1_plot_by_model_row,
       width = 5, height = 10, units = "in", dpi = 300)




##### Now create a wide version of the table 
beta1_df_wide <- jags_beta1_df %>% dplyr::select(!model_labs) %>%
             pivot_wider(
            id_cols = c(species, threshold),
            names_from = model, 
            names_prefix = "beta1_",
            values_from = c(mean, sd, `X2.5.`, `X50.`, `X97.5.`, Rhat, overlap0, direction, overlap0direction)) %>%
           # column_to_rownames(var = "species") %>%
  filter(!is.na(mean_beta1_model_S))
  

beta1_df_wide %>% filter()

mean_cols = names(beta1_df_wide)[grepl("mean_", names(beta1_df_wide))]
median_cols = names(beta1_df_wide)[grepl("X50._", names(beta1_df_wide))]
overlap0_cols = names(beta1_df_wide)[grepl("overlap0_", names(beta1_df_wide))]
non0direction_cols = names(beta1_df_wide)[grepl("overlap0direction", names(beta1_df_wide))]
sd_cols = names(beta1_df_wide)[grepl("sd", names(beta1_df_wide))]
beta1_df_wide %>% dplyr::select(all_of(overlap0_cols)) %>% rowwise() %>% mutate(count_models_overlap0 = sum(c_across(overlap0_cols)))
beta1_df_wide %>% 
  dplyr::select(all_of(overlap0_cols), species) %>% 
  mutate(count_models_overlap0 = rowSums(across(overlap0_cols)))

beta1_df_wide %>% 
  dplyr::select(all_of(non0direction_cols), species) %>% rowwise %>%
  dplyr::mutate(count_pos_models = sum(across(non0direction_cols) > 0),
         count_neg_models = sum(across(non0direction_cols) < 0),
         count_0_models = sum(across(non0direction_cols) == 0)) %>%
  dplyr::select(count_pos_models, count_neg_models, count_0_models) %>% 
  group_by(count_pos_models, count_neg_models) %>%
  tally()


get_kappa_matrix(beta1_df_wide[c(overlap0_cols)], function_name = "agreement")
get_pairs_from_matrix(as.data.frame(get_kappa_matrix(beta1_df_wide[c(overlap0_cols)], function_name = "maxwell")), var_name = "maxwell")
get_pairs_from_matrix(as.data.frame(get_kappa_matrix(beta1_df_wide[c(overlap0_cols)], function_name = "maxwell")), var_name = "maxwell")
kripp.alpha(as.matrix(beta1_df_wide[c(mean_cols)]), method = "ratio")

icc(as.matrix(beta1_df_wide[c(mean_cols)]), model = "oneway")
icc(as.matrix(beta1_df_wide[c(median_cols)]), model = "oneway")

beta1_non0direction_kappa_pairs_rank <- get_pairs_from_matrix(as.data.frame(get_kappa_matrix(beta1_df_wide[c(non0direction_cols)], 
                                                                                  function_name = "kappa")),
                                                   var_name = "kappa",
                                                   rank_fn = min_rank)

beta1_non0direction_kappa_pairs_rank %>% 
  group_by(row_name) %>% 
  summarize(mean_rank = mean(rank), mean_kappa = mean(kappa)) %>% 
  arrange(mean_rank, mean_kappa)

get_pairs_from_matrix(as.data.frame(
  get_kappa_matrix(beta1_df_wide[c(overlap0_cols)], 
                   function_name = "kappa")),
  var_name = "kappa",
  rank_fn = min_rank) %>% 
  group_by(row_name) %>% 
  summarize(mean_rank = mean(rank), mean_kappa = mean(kappa)) %>% 
  arrange(mean_rank, mean_kappa)

get_pairs_from_matrix(as.data.frame(
  get_kappa_matrix(beta1_df_wide[c(mean_cols)], 
                   function_name = "pearson_cor")),
  var_name = "pearson_cor",
  rank_fn = min_rank) %>% 
  group_by(row_name) %>% 
  summarize(mean_rank = mean(rank), mean_pearson_cor = mean(pearson_cor)) %>% 
  arrange(mean_rank, mean_pearson_cor)

get_pairs_from_matrix(as.data.frame(
  get_kappa_matrix(beta1_df_wide[c(mean_cols)], 
                   function_name = "spearman_rho")),
  var_name = "spearman_rho",
  rank_fn = min_rank) %>% 
  group_by(row_name) %>% 
  summarize(mean_rank = mean(rank), mean_spearman_rho = mean(spearman_rho)) %>% 
  arrange(mean_rank, mean_spearman_rho)


#### Multi-metric identify most discordanant model ############
beta1_ratings_list <- list(overlap0_beta1 = as.matrix(beta1_df_wide[c(overlap0_cols)]), 
                           non0direction_beta1 = as.matrix(beta1_df_wide[c(non0direction_cols)]),
                           mean_beta1 = as.matrix(beta1_df_wide[c(mean_cols)]),
                           median_beta1 = as.matrix(beta1_df_wide[c(median_cols)])
                           )
get_most_discordant_model(beta1_ratings_list, metric_func=iota)


kf_o0_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(overlap0_cols)]), 
                                      metric_func = kappam.fleiss)$metric_impact_df %>% 
  mutate(metric_function = "Kappa Fleiss", measure = "overlap0")



rob_o0_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(overlap0_cols)]), 
                                       metric_func = robinson)$metric_impact_df %>%
  mutate(metric_function = "Robinson", measure = "overlap0")


agree_o0_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(overlap0_cols)]), 
                                         metric_func = agree)$metric_impact_df %>% 
  mutate(metric_function = "Agree", measure = "overlap0", metric_with_model_removed = as.numeric(metric_with_model_removed)/100)

kf_non0dir_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(non0direction_cols)]), 
                                           metric_func = kappam.fleiss)$metric_impact_df %>% 
  mutate(metric_function = "Kappa Fleiss", measure = "Non-0 direction of beta_1")

calpha_meanbeta1_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(mean_cols)]), 
                                                   metric_func = cronbach.alpha)$metric_impact_df %>% 
  mutate(metric_function = "Chronbach alpha", measure = "Posterior mean of beta_1")

kripp_alpha_meanbeta1_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(mean_cols)]), 
                                                        metric_func = kripp.alpha)$metric_impact_df %>% 
  mutate(metric_function = "Krippendorff alpha", measure = "Posterior mean of beta_1")

icc_meanbeta1_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(mean_cols)]), 
                                              metric_func = icc)$metric_impact_df %>% 
  mutate(metric_function = "ICC", measure = "Posterior mean of beta_1")

kripp_alpha_medianbeta1_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(median_cols)]), 
                                                    metric_func = kripp.alpha)$metric_impact_df %>% 
  mutate(metric_function = "Krippendorff alpha", measure = "Posterior median of beta_1")

icc_medianbeta1_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(median_cols)]), 
                                                metric_func = icc)$metric_impact_df %>% 
  mutate(metric_function = "ICC", measure = "Posterior median of beta_1")

calpha_medianbeta1_df <- get_most_discordant_model(as.matrix(beta1_df_wide[c(median_cols)]), 
                                                metric_func = cronbach.alpha)$metric_impact_df %>% 
  mutate(metric_function = "Chronbach alpha", measure = "Posterior median of beta_1")

cronbach.alpha(beta1_df_wide[c(median_cols)], standardized = TRUE, CI = TRUE, 
               probs = c(0.025, 0.975), B = 1000, na.rm = FALSE)


icc_sdbeta1_df  <- get_most_discordant_model(as.matrix(beta1_df_wide[c(sd_cols)]), 
                                                                  metric_func = icc)$metric_impact_df %>% 
  mutate(metric_function = "ICC", measure = "Posterior sd of beta_1")

combined_discordant_df <- rbind.data.frame(kf_o0_df,rob_o0_df, agree_o0_df, kf_non0dir_df,
                                           icc_meanbeta1_df, icc_medianbeta1_df, icc_sdbeta1_df)

cronbach.alpha(beta1_df_wide[c(median_cols)], standardized = TRUE, CI = TRUE, 
               probs = c(0.025, 0.975), B = 1000, na.rm = FALSE)

combined_discordant_df$combined_metric <- paste(combined_discordant_df$metric_function, combined_discordant_df$measure, sep = ": \n")


# Add a count of duplicates to the data
df_count <- combined_discordant_df %>%
  group_by(combined_metric, concordance_rank) %>%
  mutate(n = n()) %>%
  ungroup()

# Plot using conditional logic
ggplot(df_count, aes(x, y)) +
  # Plot points with no ties using geom_point
  geom_point(data = filter(df_count, n == 1)) +
  # Plot jittered points for ties, using the position argument
  geom_point(data = filter(df_count, n > 1),
             position = position_jitter(width = 0.1, height = 0.1))


my_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
               "#e6ab02", "#a6761d", "#666666", "#e78ac3", "#984ea3")

rank_plot_beta1 <- ggplot(df_count, aes(x = combined_metric, y = concordance_rank, group = model_labels, color = model_labels)) +
  geom_line(size = .3) +
  geom_point(data = filter(df_count, n == 1), size = 6) +
  geom_jitter(data = filter(df_count, n > 1), size = 6, width=0.1, height = 0.1) +
  geom_text(aes(label = concordance_rank), vjust = 0.6, size = 3.5, color = 'white') +
  scale_color_manual(values = my_colors) +
  scale_y_reverse(breaks = 1:nrow(combined_discordant_df)) +
  labs(x = "", y = "") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                legend.title = element_blank(), 
        legend.position = 'bottom',
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )  + 
  # scale_y_reverse(
  #   breaks = 1:n_distinct(combined_discordant_df$model_labels),
  #   labels = combined_discordant_df %>% filter(combined_metric == 'Posterior median of beta_1') %>% arrange(concordance_rank) %>% pull(model_labels)
  #   # sec.axis = sec_axis(
  #   #   trans = I,
  #   #   breaks = 1:n_distinct(df$model_labels),
  #   #   labels = df %>% filter(year == 'year4') %>% 
  #   #     arrange(rank) %>% pull(car))
  #   # 
  #   ) +  
  scale_x_discrete(expand = c(0,.1))

ggsave("discordant_metrics_rank_beta1.png", plot = rank_plot_beta1, width = 10, height = 5, units = "in", dpi = 300)


library(tidytext)
library(ggh4x)

df_count <- df_count %>%
  mutate(model_labels_ordered = reorder_within(model_labels, as.numeric(metric_with_model_removed), combined_metric, fun = "identity"))

ggplot(df_count, aes(x = model_labels_ordered, 
                     y = as.numeric(metric_with_model_removed),
                     fill = model_labels),
       ) +
  scale_x_reordered() +
  geom_bar(stat="identity") + 
  facet_nested(~ measure + metric_function, scales = "free") +  theme_bw()+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #   legend.title = element_blank(), legend.position = 'bottom',    panel.border = element_blank()
  # ) +
  labs(x = "Models ordered from least to most discordant", y = "Metric value\nwith model removed")

metric_fn_labeller <- function(variable, value) {
  if (variable == "combined_metric") {
    return(switch(value,
           "Kappa Fleiss: \nNon-0 direction of beta_1" = expression(paste0("Kappa Fleiss: \nSign of non-0 ", beta[1])),
           "ICC: \nPosterior mean of beta_1" = expression(paste0("ICC: \nMean ", beta[1])),
           "ICC: \nPosterior sd of beta_1" = expression(paste0("ICC: \nSD", beta[1])),
           "ICC: \nPosterior median of beta_1" = expression(paste0("ICC: \nMedian ", beta[1])),
           "Kappa Fleiss: \noverlap0" = expression(paste0("Kappa Fleiss: \n95% HDI of \n", beta[1], " overlaps 0")),
           "Robinson: \noverlap0" = expression(paste0("Robinson: \n95% HDI of \n", beta[1], " overlaps 0")),
           "Agree: \noverlap0" = expression(paste0("Agree: \n95% HDI of \n", beta[1], " overlaps 0"))
           ))
  } else{
      return(value)
    }
}

ggplot(df_custom, aes(x, y)) +
  geom_point() +
  facet_wrap(~ type, labeller = labeller(my_labeller, default = label_parsed))


combined_labels <- c(
  "Kappa Fleiss: \nNon-0 direction of beta_1" = expression(paste0("Kappa Fleiss: \nSign of non-0 ", beta[1])),
  "ICC: \nPosterior mean of beta_1" = expression(paste0("ICC: \nMean ", beta[1])),
  "ICC: \nPosterior sd of beta_1" = expression(paste0("ICC: \nSD", beta[1])),
  "ICC: \nPosterior median of beta_1" = expression(paste0("ICC: \nMedian ", beta[1])),
  "Kappa Fleiss: \noverlap0" = expression(paste0("Kappa Fleiss: \n95% HDI of \n", beta[1], " overlaps 0")),
  "Robinson: \noverlap0" = expression(paste0("Robinson: \n95% HDI of \n", beta[1], " overlaps 0")),
  "Agree: \noverlap0" = expression(paste0("Agree: \n95% HDI of \n", beta[1], " overlaps 0"))
)
     
discordant_metrics_plot <- ggplot(df_count, aes(x = model_labels_ordered, 
                     y = as.numeric(metric_with_model_removed),
                     fill = model_labels)) +
  scale_x_reordered() +
  geom_bar(stat="identity") + 
  geom_text(aes(label = round(as.numeric(metric_with_model_removed), 2)), vjust = 1.5, size=3) +
  facet_wrap(~ as.factor(combined_metric), ncol = 1, scales = "free", strip.position = "top") +  
  theme_bw()+
  theme(strip.text.y = element_text(angle = 0)) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #   legend.title = element_blank(), legend.position = 'bottom',    panel.border = element_blank()
  # ) +
  labs(x = "Models ordered from least to most discordant",
       y = "Metric value\nwith model removed",
       fill = "Model")

ggsave("discordant_metrics_beta1.png", plot = discordant_metrics_plot, width = 5, height = 10, units = "in", dpi = 300)

all_b1_col_names = names(beta1_df_wide)
rhat_col_names <- all_b1_col_names[grepl(paste0("Rhat_beta1_"), all_b1_col_names)] 
beta1_df_wide %>% filter(if_any(rhat_col_names, ~ .x > 1.5) )
                        
                         
                         
                         
                         
filter_non_converged_beta_results <- function(results_df, 
                                         kos = c("model_A", "model_AS", 
                                                 "model_H", "model_HA", 
                                                 "model_HAS", "model_HS", 
                                                 "model_S"), 
                                         rhat_tol = 1.5) {
  converged_median_cols <- c()
  converged_mean_cols <- c()
  for (ko in kos){
    all_col_names = names(results_df)
    rhat_col_name <- all_col_names[grepl(paste0("Rhat_beta1_",ko,"$"), all_col_names)]
    median_col_name <- all_col_names[grepl(paste0("X50._beta1_",ko,"$"), all_col_names)]
    mean_col_name <- all_col_names[grepl(paste0("mean_beta1_",ko,"$"), all_col_names)]
    if(sum(occ_df[[rhat_col_name]] > rhat_tol, na.rm=TRUE) <= n_sites_tolerance) {
      converged_median_cols <- c(converged_median_cols, median_col_name)
      converged_mean_cols <- c(converged_mean_cols, mean_col_name)
    }
  }
  return(list(median_cols = converged_median_cols, mean_cols = converged_mean_cols))
}

lbrary(lme4)
lmer(mean ~ 1 + model + (1 | species), data = jags_beta1_df)
mod1 <- lmer(mean ~ model + (1 | species), data = jags_beta1_df)
ranef(mod1)

mod2 <- glmer(overlap0 ~ model + (1 | species), data = jags_beta1_df, family = "binomial")
ranef(mod2)

beta1_df_wide[c(median_cols)]

library(corrplot)

beta_agree_matrix <- get_kappa_matrix(beta1_df_wide[c(overlap0_cols)], 
                 function_name = "agreement")
corrplot(beta_agree_matrix,addCoef.col = 'grey85', order = 'AOE')

beta_agree_matrix <- get_kappa_matrix(beta1_df_wide[c(overlap0_cols)], 
                                      function_name = "agreement")
corrplot(cor(beta1_df_wide[c(mean_cols)]),addCoef.col = 'grey85', order = 'AOE')

pairs(beta1_df_wide[c(mean_cols)])

library(GGally)
add_y_x_line <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(...) + # Include original plot elements if desired
    geom_abline(intercept = 0, slope =  1, linetype = "dashed", color = "blue")  + 
    scale_x_continuous(limits = c(-4.5, 4)) + 
    scale_y_continuous(limits = c(-4.5, 4))
  return(p)
}

model_labeller <- function(label_cols) {
  # Add a prefix
  str_split_i(label_cols, "_", i = -1)
  # Wrap long labels (requires 'label_wrap_gen' from ggplot2)
}

mean_cols %>% str_split_i("_", -1)

ggpairs(beta1_df_wide, columns = mean_cols, lower = list(continuous = add_y_x_line), 
        labeller=as_labeller(model_labeller))
      
min(beta1_df_wide[c(mean_cols)])
max(beta1_df_wide[c(mean_cols)])
        
