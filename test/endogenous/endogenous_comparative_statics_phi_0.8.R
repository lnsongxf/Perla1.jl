library(ggplot2)
library(cowplot)
library(dplyr)

# Plot setup
SAVE.INDIVIDUAL.PLOTS <- FALSE
TEXT.SIZE <- 15
PLOT.SIZE.WIDTH <- 5
PLOT.SIZE.HEIGHT <- 5
PLOT.NAME <- "endogenous_comparative_statics_phi_0.8"
theme_set(theme_bw())

# plot generation definitions
GenerateSharePlot <- function (param.name, df) {
  df.plot <- df %>% filter(change_by == param.name)
  
  # generate plot
  scale.factor <- 25
  plot.variable <- ggplot(df.plot, aes_string(x=param.name)) +
    geom_line(aes(y=mu, col="mu")) +
    geom_line(aes(y=theta*scale.factor, col="theta"))  +
    scale_y_continuous(name="mu", sec.axis=sec_axis(~./scale.factor, name="theta")) +
    scale_color_manual(values=c("#3175BC", "#E34F28")) +
    labs(colour="Variable") +
    guides(size = FALSE) +
    theme(text = element_text(size=TEXT.SIZE))
  
  # remove legend if needed
  if (REMOVE.LEGEND)
    plot.variable <- plot.variable + theme(legend.position = "none")
  
  # save plot
  plot.variable.name.png <- paste0(PLOT.NAME, "-theta-and-mu-by-", param.name, ".png")
  plot.variable
  if (SAVE.INDIVIDUAL.PLOTS)
    ggsave(plot.variable.name.png, width = PLOT.SIZE.WIDTH, height = PLOT.SIZE.HEIGHT)
  
  return (plot.variable)
}

GenerateVariablePlot <- function(param.name, df) {
  df.plot <- df %>% filter(change_by == param.name)
  
  # generate plot
  scale.factor <- max(df.plot$profit_share) / max(df.plot$Q_by_B)
  plot.variable <- ggplot(df.plot, aes_string(x=param.name)) +
    geom_line(aes(y=profit_share, col="Profit share")) +
    geom_line(aes(y=Q_by_B*scale.factor, col="Q/B"))  +
    scale_y_continuous(name="Profit share", sec.axis=sec_axis(~./scale.factor, name="Q/B")) + 
    scale_color_manual(values=c("#3175BC", "#E34F28")) +
    labs(colour="Variable") +
    guides(size = FALSE) +
    theme(text = element_text(size=TEXT.SIZE))
  
  # remove legend if needed
  if (REMOVE.LEGEND)
    plot.variable <- plot.variable + theme(legend.position = "none")
  
  # save plot
  plot.variable.name.png <- paste0(PLOT.NAME, "-variables-by-", param.name, ".png")
  plot.variable
  if (SAVE.INDIVIDUAL.PLOTS)
    ggsave(plot.variable.name.png, width = PLOT.SIZE.WIDTH, height = PLOT.SIZE.HEIGHT)
  
  return (plot.variable)
}

# load df
df <- read.csv(paste0(PLOT.NAME, ".csv"))
df <- df %>% 
  mutate(Q_by_B = Q / B)

# generate plots
REMOVE.LEGEND <- TRUE 
plots.share <- lapply(levels(df$change_by), GenerateSharePlot, df = df)
plots.variable <- lapply(levels(df$change_by), GenerateVariablePlot, df = df)

REMOVE.LEGEND <- FALSE
plot.share <- GenerateSharePlot(levels(df$change_by)[1], df = df)
plot.variable <-  GenerateVariablePlot(levels(df$change_by)[1], df = df)

# plots together
plot_grid(plotlist=c(c(plots.share, list(get_legend(plot.share + theme(legend.justification ="left")))), 
                     c(plots.variable, list(get_legend(plot.variable+ theme(legend.justification ="left"))))), ncol=(length(levels(df$change_by)) + 1)) +
  theme(plot.margin = margin(0, -PLOT.SIZE.WIDTH*0.7, 0, 0, "in")) # remove margin on the legend
ggsave(paste0(PLOT.NAME, ".png"), width = (PLOT.SIZE.WIDTH*4), height = (PLOT.SIZE.HEIGHT*2))
