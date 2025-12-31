create_plot <- function(data, chart_name, UCL, stat_name) {
  
  shifts <- 1:nrow(data)
  ARL <- data[[chart_name]]
  
  plot_data <- data.frame(Shifts = shifts, ARL = ARL)
  
  ggplot(plot_data, aes(x = Shifts, y = ARL)) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 2.2) +
    
    # --- UCL line ---
    geom_hline(yintercept = UCL, linewidth = 0.8) +
    geom_text(aes(x = max(Shifts) + 0.01, y = UCL,
                  label = "UCL"),
              size = 5, hjust = 0, vjust = -0.5) +
    
    # --- Labels ---
    labs(x = "Sample Number", y = stat_name) +
    
    # --- Auto axis scaling ---
    scale_x_continuous(breaks = seq(0, max(shifts), by = 2),
                       expand = expansion(mult = c(0.02, 0.07))) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.10))  # <-- auto space above/below
    ) +
    
    # --- CLEAN THEME (NO BACKGROUND) ---
    theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      
      legend.position = "none",
      axis.ticks.length = unit(0.30, "cm")
    )
}

p1 <- create_plot(plot_df, "E_VWAB", UCL = 2.20, stat_name = bquote(italic(EL["V,A"])))
print(p1)

p2 <- create_plot(plot_df, "E_VWM", UCL = 2.195, stat_name = bquote(italic(EL["V,M"])))
print(p2)

p3 <- create_plot(plot_df, "E_WRSAB", UCL = 2.205, stat_name = bquote(italic(EL["W,A"])))
print(p3)

p4 <- create_plot(plot_df, "E_C", UCL = 1.10, stat_name = bquote(italic(EC)))
print(p4)
