
#################################################
### VISUALISATION FUNCTIONS FOR SEIR EPIDEMICS ##
#################################################

library(gridExtra)

#~~~~~~~~~~~~~#
# TRACE PLOTS #
#~~~~~~~~~~~~~#

gg_trace_plot <- function(results, params_true, burn_in, annotate_xy, x_limits){
  
  df = data.frame(samples = results[, 4], b = results[, 1], g = results[, 2], d = results[, 3])
  
  num_samples = nrow(df)
  
  df_burn_in <- df %>% 
    mutate(burnin = c(rep("yes", burn_in), rep("no", (num_samples - burn_in))))
  
  beta_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = b)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    annotate("text", y = annotate_xy[2], x = annotate_xy[1], ##### HERE IS CHANGE VALUE
             parse = TRUE, size = 5,
             label = as.expression(bquote(beta~"="~.(params_true[[1]])))) + ##### HERE IS CHANGE SYMBOL AND VALUE
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(x_limits[1], x_limits[2])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ beta)) + 
    guides(colour = "none") 
  
  
  
  g_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = g)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[2]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    annotate("text", y = annotate_xy[4], x = annotate_xy[3], ##### HERE IS CHANGE VALUE
             parse = TRUE, size = 5,
             label = as.expression(bquote(delta~"="~.(params_true[[2]])))) + ##### HERE IS CHANGE SYMBOL AND VALUE
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(x_limits[3], x_limits[4])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ delta)) + 
    guides(colour = "none") 
  
  
  d_plot = df_burn_in %>% 
    ggplot(aes(x = samples, y = d)) +
    geom_line(size  = 0.2, aes(colour = burnin)) + 
    scale_colour_manual(values = c("dodgerblue", "#fd8f24")) +
    geom_hline(aes(yintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = '#4f5157') +
    annotate("text", y = annotate_xy[6], x = annotate_xy[5], ##### HERE IS CHANGE VALUE
             parse = TRUE, size = 5,
             label = as.expression(bquote(gamma~"="~.(params_true[[3]])))) + ##### HERE IS CHANGE SYMBOL AND VALUE
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(x_limits[5], x_limits[6])) +
    labs(
      x = "Sample index", 
      y = "Value", 
      title = expression("Traceplot of parameter samples" ~ gamma)) + 
    guides(colour = "none") 
  
  trace_plots = grid.arrange(beta_plot, g_plot, d_plot, nrow = 3)
  
  # use ggsave() to save it, e.g.
  # ggsave(filename = "F_augs_trace.png", width = 10, height = 8) ##### HERE IS CHANGE NAME
  
  return(trace_plots)
}


#~~~~~~~~~~~~#
# HISTOGRAMS #
#~~~~~~~~~~~~#

gg_hist_plot <- function(results, params_true, burn_in, annotate_xy, x_limits){
  
  df = data.frame(samples = results[, 4], b = results[, 1], g = results[, 2], d = results[, 3])
  
  num_samples = nrow(df)
  
  df_burn_in <- df %>% 
    mutate(burnin = c(rep("yes", burn_in), rep("no", (num_samples - burn_in))))
  
  beta_plot = df_burn_in %>% 
    ggplot(aes(x = b, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[1]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(x_limits[1], x_limits[2])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    annotate("text", y = annotate_xy[2], x = annotate_xy[1], ##### HERE IS CHANGE VALUE
             parse = TRUE, size = 5,
             label = as.expression(bquote(beta~"="~.(params_true[[1]])))) + ##### HERE IS CHANGE SYMBOL AND VALUE
    labs(
      x = expression("Samples of parameter" ~ beta), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  g_plot = df_burn_in %>% 
    ggplot(aes(x = g, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[2]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(x_limits[3], x_limits[4])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    annotate("text", y = annotate_xy[4], x = annotate_xy[3], ##### HERE IS CHANGE VALUE
             parse = TRUE, size = 5,
             label = as.expression(bquote(delta~"="~.(params_true[[2]])))) + ##### HERE IS CHANGE SYMBOL AND VALUE
    labs(
      x = expression("Samples of parameter" ~ delta), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  d_plot = df_burn_in %>% 
    ggplot(aes(x = d, y = ..density..)) +
    geom_histogram(fill = "#f5c04a", colour = "grey15", alpha = 0.85) +
    geom_vline(aes(xintercept = params_true[3]), size = 1.1, linetype = 2, 
               colour = '#4f5157') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), 
                       limits = c(x_limits[5], x_limits[6])) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)
    ) +
    annotate("text", y = annotate_xy[6], x = annotate_xy[5], ##### HERE IS CHANGE VALUE
             parse = TRUE, size = 5,
             label = as.expression(bquote(gamma~"="~.(params_true[[3]])))) + ##### HERE IS CHANGE SYMBOL AND VALUE
    labs(
      x = expression("Samples of parameter" ~ gamma), ##### HERE IS CHANGE SYMBOL
      y = "Density") 
  
  
  hist_plots = grid.arrange(beta_plot, g_plot, d_plot, nrow = 1)
  
  return(hist_plots)
}





gg_contour_plot <- function(results, params_true, burn_in, b_limits, g_limits, d_limits){
  
  df_burn_in = data.frame(b = results[-(1:burn_in), 1], g = results[-(1:burn_in), 2], d = results[-(1:burn_in), 3])
  
  
  
  
  hexplot_init = ggplot(df_burn_in, aes(x=b, y=g) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50) +
    scale_fill_continuous(type = "viridis")
  
  meta_data = ggplot_build(hexplot_init)$data
  
  mean_x = meta_data[[2]][which.max(meta_data[[2]]$density),]$x
  mean_y = meta_data[[2]][which.max(meta_data[[2]]$density),]$y
  mean_beta = signif(mean_x, digits = 3)
  mean_gamma = signif(mean_y, digits = 3) 
  
  
  rm(hexplot_init)
  
  
  
  beta_plot = df_burn_in %>% 
    ggplot(aes(x = b, y = g)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[2]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$b), max(df_burn_in$b))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$g), max(df_burn_in$g))) +
    coord_cartesian(xlim = c(b_limits[1], b_limits[2]),  ylim = c(g_limits[1], g_limits[2])) +
    labs(
      x = as.expression(bquote(beta)), 
      y = as.expression(bquote(delta)), 
      title = expression("Contour plot of parameter samples" ~ beta ~ "vs." ~ delta)) + 
    guides(colour = "none") 
  
  
  
  
  
  
  
  
  hexplot_init_1 = ggplot(df_burn_in, aes(x=b, y=d) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50)
  
  meta_data_1 = ggplot_build(hexplot_init_1)$data
  
  mean_x_1 = meta_data_1[[2]][which.max(meta_data_1[[2]]$density),]$x
  mean_y_1 = meta_data_1[[2]][which.max(meta_data_1[[2]]$density),]$y
  mean_beta_1 = signif(mean_x_1, digits = 3)
  mean_delta_1 = signif(mean_y_1, digits = 3) 
  
  
  g_plot = df_burn_in %>% 
    ggplot(aes(x = b, y = d)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x_1), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y_1), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[1]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$b), max(df_burn_in$b))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$d), max(df_burn_in$d))) +
    coord_cartesian(xlim = c(b_limits[1], b_limits[2]),  ylim = c(d_limits[1], d_limits[2])) +
    labs(
      x = as.expression(bquote(beta)), 
      y = as.expression(bquote(gamma)), 
      title = expression("Contour plot of parameter samples" ~ beta ~ "vs." ~ gamma)) + 
    guides(colour = "none") 
  
  
  
  
  
  
  
  hexplot_init_2 = ggplot(df_burn_in, aes(x=g, y=d) ) +
    geom_point(color = NA, fill = NA) +
    geom_hex(bins = 50)
  
  meta_data_2 = ggplot_build(hexplot_init_2)$data
  
  mean_x_2 = meta_data_2[[2]][which.max(meta_data_2[[2]]$density),]$x
  mean_y_2 = meta_data_2[[2]][which.max(meta_data_2[[2]]$density),]$y
  mean_gamma_2 = signif(mean_x_2, digits = 3)
  mean_delta_2 = signif(mean_y_2, digits = 3) 
  
  
  
  
  
  
  
  d_plot = df_burn_in %>% 
    ggplot(aes(x = g, y = d)) +
    stat_density2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_continuous(type = "viridis") +
    geom_vline(aes(xintercept = mean_x_2), size = 0.7, linetype = 2, 
               colour = 'yellow')+
    geom_hline(aes(yintercept = mean_y_2), size = 0.7, linetype = 2, 
               colour = 'yellow') +  
    geom_vline(aes(xintercept = params_true[2]), size = 0.7, linetype = 2, 
               colour = 'red')+
    geom_hline(aes(yintercept = params_true[3]), size = 0.7, linetype = 2, 
               colour = 'red') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(min(df_burn_in$g), max(df_burn_in$g))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(min(df_burn_in$d), max(df_burn_in$d))) +
    coord_cartesian(xlim = c(g_limits[1], g_limits[2]),  ylim = c(d_limits[1], d_limits[2])) +
    labs(
      x = as.expression(bquote(delta)), 
      y = as.expression(bquote(gamma)), 
      title = expression("Contour plot of parameter samples" ~ delta ~ "vs." ~ gamma)) + 
    guides(colour = "none") 
  
  trace_plots = grid.arrange(beta_plot, g_plot, d_plot, nrow = 3)
  
  # use ggsave() to save it, e.g.
  # ggsave(filename = "F_augs_trace.png", width = 10, height = 8) ##### HERE IS CHANGE NAME
  
  
  return(trace_plots)
}











