library(tidyverse)

docker <- read_csv('./results-docker.csv') %>% mutate(environment = 'docker')
singularity <- read_csv('./results-singularity.csv') %>% mutate(environment = 'singularity')
native <- read_csv('./results-native.csv') %>% mutate(environment = 'native')

results <- bind_rows(docker, singularity, native);
results <- results %>% 
  group_by(environment, size) %>%
  summarize(
    samples = n(),
    average = mean(time),
    stdDeviation = sd(time),
    stdError = 3*stdDeviation/sqrt(samples))

default_theme <- function() {
  ret <- list();
  ret[[length(ret)+1]] <- theme (
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.spacing = unit(1, "mm"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.spacing = unit(0, "pt"),
    legend.box.margin = margin(0,0,0,0),
    legend.title = element_blank());
  return(ret);
}

g <- ggplot(results,aes(x=size, y=average)) +
  geom_line(aes(col = environment), alpha = 0.2) +
  geom_point(aes(col = environment), size = 3) +
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, color=environment, group=environment), width = 0.3) +
  theme_bw(base_size=12) +
  scale_y_continuous(trans='log2') + 
  #ylim(0,NA) +
  scale_x_continuous(trans="log2") + 
  ylab('Average latency (ms)') +
  xlab('Message size (bytes)') +
  scale_color_grey() +
  default_theme()

plot(g)
