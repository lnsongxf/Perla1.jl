library(ggplot2)
library(dplyr)
library(cowplot)

# Plot setup
SAVE.INDIVIDUAL.PLOTS <- FALSE
TEXT.SIZE <- 15
PLOT.SIZE.WIDTH <- 5
PLOT.SIZE.HEIGHT <- 5
theme_set(theme_bw())

# plot generation definitions
GenerateSharePlot <- function (param.name, df) {
  df.plot <- df %>% filter(change_by == param.name)
  
  # generate plot
  plot.share <- ggplot(df.plot, aes_string(x=param.name)) +
    geom_area(aes(y=capital_share+profit_share, fill = "Profit share")) +
    geom_area(aes(y=capital_share, fill = "Capital share")) + 
    scale_fill_manual(values=c("#3175BC", "#E34F28")) +
    labs(y="Shares", fill="Share") +
    ylim(0.0, 0.52) +
    theme(text = element_text(size=TEXT.SIZE))
  
  # remove legend if needed
  if (REMOVE.LEGEND)
    plot.share <- plot.share + theme(legend.position = "none")
  
  # save plot
  plot.share.name.png <- paste0("shares-by-", param.name, ".png")
  plot.share
  if (SAVE.INDIVIDUAL.PLOTS)
    ggsave(plot.share.name.png, width = PLOT.SIZE.WIDTH, height = PLOT.SIZE.HEIGHT)
  
  return (plot.share)
}

GenerateVariablePlot <- function(param.name, df) {
  df.plot <- df %>% filter(change_by == param.name)
  
  # generate plot
  scale.factor <- max(df.plot$Q_by_B) / max(df.plot$Z)
  plot.variable <- ggplot(df.plot, aes_string(x=param.name)) +
    geom_line(aes(y=Q_by_B, col="Q/B")) +
    geom_line(aes(y=Z*scale.factor, col="Z"))  +
    scale_y_continuous(name="Q/B", sec.axis=sec_axis(~./scale.factor, name="Z")) + 
    scale_color_manual(values=c("#3175BC", "#E34F28")) +
    labs(colour="Variable") +
    guides(size = FALSE) +
    theme(text = element_text(size=TEXT.SIZE))
  
  # remove legend if needed
  if (REMOVE.LEGEND)
    plot.variable <- plot.variable + theme(legend.position = "none")
  
  # save plot
  plot.variable.name.png <- paste0("variables-by-", param.name, ".png")
  plot.variable
  if (SAVE.INDIVIDUAL.PLOTS)
    ggsave(plot.variable.name.png, width = PLOT.SIZE.WIDTH, height = PLOT.SIZE.HEIGHT)
  
  return (plot.variable)
}

# load df
df <- read.csv("exogenous_comparative_statics.csv")
df <- df %>% 
  mutate(Q_by_B = Q / B)

# rescale Z by baseline and drop baseline 
Z.baseline <- (df %>% filter(change_by == "baseline"))$Z
df <- df %>%
  mutate(Z = Z / Z.baseline) %>%
  filter(change_by != "baseline") %>% 
  droplevels()

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
ggsave("exogenous_comparative_statics.png", width = (PLOT.SIZE.WIDTH*4), height = (PLOT.SIZE.HEIGHT*2))
