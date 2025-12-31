library(DescTools)
library(ggplot2)
source("Distance_Statistics.R")
df <- read.csv("Piston Ring Diameter Data.csv")
head(df)

Phase_1 = df[df$Phase=="I",paste0("X", 1:5)]
Phase_2 = df[df$Phase=="II",paste0("X", 1:5)]
u = unlist(lapply(paste0("X", 1:5), function(x) Phase_1[[x]]))
print(paste0("m=",length(u)))


# Define UDF to process and update plot_df
generate_plot_df <- function(phase_2, u,lambda=0.05) {
  # Initialize plot_df with default values
  plot_df <- data.frame(
    t = 0,
    ZWRS = 0,
    ZVW = 0,
    ZAB = 0,
    ZM = 0,
    SC = 0,
    WRSAB = 0,
    WRSM = 0,
    MWAB = 0,
    MWM = 0,
    VWAB = 0,
    VWM = 0,
    PvalMW = 0,
    PvalVW = 0,
    PvalAB = 0,
    PvalM = 0,
    PvalMN = 0
  )
  
  # Populate plot_df with values calculated for each sample
  for (t in 1:nrow(phase_2)) {
    v <- c(as.vector(unlist(phase_2[t,]))) #c(t,t-1)
    print(length(v))
    
    # Update plot_df with computed statistics for sample t
    plot_df[t, ] <- c(
      t,
      ZWRS(u,v),
      ZVW(u,v),
      ZAB(u,v),
      ZM(u,v),
      SC(u, v),
      WRSAB(u, v),
      WRSM(u, v),
      MWAB(u, v),
      MWM(u, v),
      VWAB(u, v),
      VWM(u, v),
      wilcox.test(u, v)$p.value,
      VanWaerdenTest(list(v, u))$p.value,
      ansari.test(u, v)$p.value,
      mood.test(u, v)$p.value,
      mood.test(u - median(u), v - median(v))$p.value
    )
  }
  
  for( t in 1:nrow(Phase_2)){
    if(t==1){
      plot_df[t,"E_C"] = max(1,lambda*plot_df[t,"SC"]+(1-lambda)*1)
      plot_df[t,"E_WRSAB"] = max(2,lambda*plot_df[t,"WRSAB"]+(1-lambda)*2)
      plot_df[t,"E_WRSM"] = max(2,lambda*plot_df[t,"WRSM"]+(1-lambda)*2)
      plot_df[t,"E_MWAB"] = max(2,lambda*plot_df[t,"MWAB"]+(1-lambda)*2)
      plot_df[t,"E_MWM"] = max(2,lambda*plot_df[t,"MWM"]+(1-lambda)*2)
      plot_df[t,"E_VWAB"] = max(2,lambda*plot_df[t,"VWAB"]+(1-lambda)*2)
      plot_df[t,"E_VWM"] = max(2,lambda*plot_df[t,"VWM"]+(1-lambda)*2)
    }else{
      plot_df[t,"E_C"] = max(1,lambda*plot_df[t,"SC"]+(1-lambda)*plot_df[t-1,"E_C"])
      plot_df[t,"E_WRSAB"] = max(2,lambda*plot_df[t,"WRSAB"]+(1-lambda)*plot_df[t-1,"E_WRSAB"])
      plot_df[t,"E_WRSM"] = max(2,lambda*plot_df[t,"WRSM"]+(1-lambda)*plot_df[t-1,"E_WRSM"])
      plot_df[t,"E_MWAB"] = max(2,lambda*plot_df[t,"MWAB"]+(1-lambda)*plot_df[t-1,"E_MWAB"])
      plot_df[t,"E_MWM"] = max(2,lambda*plot_df[t,"MWM"]+(1-lambda)*plot_df[t-1,"E_MWM"])
      plot_df[t,"E_VWAB"] = max(2,lambda*plot_df[t,"VWAB"]+(1-lambda)*plot_df[t-1,"E_VWAB"] )
      plot_df[t,"E_VWM"] = max(2,lambda*plot_df[t,"VWM"]+(1-lambda)*plot_df[t-1,"E_VWM"])
    }
  }
  
  return(plot_df)
}

plot_df <- generate_plot_df(Phase_2, u,lambda=0.01)
print(plot_df)


# Define the plot function (UDF)
create_plot <- function(data, chart_name,lower_limit, upper_limit,stat_name) {
  dataset_name <- ""
  shifts <- 1:nrow(data)
  ARL <- data[[chart_name]]
  chart <- rep("Statistic", length(shifts))
  plot_data <- data.frame(chart, shifts, ARL)
  colnames(plot_data) <- c("Charts", "Shifts", "ARL")
  
  # Generate the plot
  plot <- ggplot(plot_data, aes(x = Shifts, y = ARL)) +
    geom_line(color = "black", size = 1.1) +
    geom_point(color = "black", size = 2.2) +
    theme_bw() +
    geom_hline(yintercept = upper_limit, size = 0.9) +
    geom_hline(yintercept = lower_limit) +
    theme(
      legend.position = c(0.2, 0.9),
      legend.title = element_blank(),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      text = element_text(size = 20),
      axis.text = element_text(size = 10),
      axis.ticks.length = unit(0.45, "cm")
    ) +
    xlab("Sample Number") +
    ylab(stat_name) +
    scale_x_continuous(breaks = seq(0, nrow(data), by = 2)) +
    scale_y_continuous(limits = c(floor(min(ARL)), ceiling(max(ARL))),
                       breaks = seq(floor(min(ARL)), ceiling(max(ARL)), by = 0.5)) +
    geom_text(aes(x = nrow(data) + 1, label = "UCL", y = upper_limit - 0.8), size = 5)
  return(plot)
}


# Generate and display the plot
ARL_graph <- create_plot(plot_df,"E_VWAB",0,2.200,bquote(italic(EL["V,A"])))
print(ARL_graph)

ARL_graph <- create_plot(plot_df,"E_VWM",0,2.195,bquote(italic(EL["V,M"])))
print(ARL_graph)

ARL_graph <- create_plot(plot_df,"E_WRSAB",0,2.205,bquote(italic(EL["W,A"])))
print(ARL_graph)

ARL_graph <- create_plot(plot_df,"E_C",0,1.100,bquote(italic('EC')))
print(ARL_graph)


# Save the final data set
output_file <- "Latest Example Piston Ring Diameter Statistics for n=5 .csv"
write.csv(plot_df, output_file)
