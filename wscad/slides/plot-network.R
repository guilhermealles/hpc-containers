library(tidyverse)
library(RColorBrewer)
# Read & mutate data
df.pingpong <- read_csv('./results/ping-pong/results.csv', col_types=cols(
  size = col_integer(),
  time = col_double(),
  environment = col_character())) %>% 
  group_by(environment, size) %>%
  summarize(
    samples = n(),
    average = mean(time),
    stdDeviation = sd(time),
    stdError = 3*stdDeviation/sqrt(samples)) %>%
  ungroup() %>%
  mutate(environment = factor(environment,
                              levels=c("native", "singularity", "docker"),
                              labels=c("Native", "Singularity", "Docker")))

# Compute color scale
colors <- brewer.pal(3,"Set1")
myColors <- c(colors[1], colors[3], colors[2])
names(myColors) <- levels(df.pingpong$environment)
colScale <- scale_colour_manual(name = "environment",values = myColors)

# Plot
df.pingpong %>%
  ggplot(aes(x=size, y=average)) +
  geom_line(aes(col = environment), alpha = 0.2) +
  geom_point(aes(col = environment), size = 1) +
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, col=environment, group=environment), width = 0.2) +
  theme_bw(base_size=12) +
  colScale +
  scale_y_log10(breaks=c(0.1, 1, 2, 4, 8, 16, 32)) +
  scale_x_log10(breaks=2^seq(0,20)) +
  ylab('Average latency (ms in logscale)') +
  xlab('Message size (bytes)') +
  theme_bw(base_size = 13) +
  theme (plot.margin = unit(c(0,0,0,0), "cm"),
         legend.spacing = unit(x = c(0, 0, 0, 0), units = 'mm'),
         legend.position = "top",
         legend.justification = "left",
         legend.box.spacing = unit(0, "pt"),
         legend.box.margin = margin(0,0,0,0),
         legend.title = element_blank(),
         axis.text.x = element_text(angle=55, hjust=1))
