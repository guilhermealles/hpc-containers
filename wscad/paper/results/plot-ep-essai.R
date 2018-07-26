library(tidyverse);

results_ep <- read_csv('./nas/results.csv') %>% mutate(experiment = 'NAS EP')
results_ondes3d <- read_csv('./ondes3d/results.csv') %>% mutate(experiment = 'Ondes3D')

results <-bind_rows(results_ep, results_ondes3d)
results <- results %>%
  mutate(time=time/1000) %>%
  group_by(experiment, environment, parallelism) %>%
  summarize(
    samples = n(),
    average = mean(time),
    stdDeviation = sd(time),
    stdError = 3*stdDeviation/sqrt(samples)
  );

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

g <- ggplot(results, aes(x = parallelism, y = average)) + 
  scale_x_continuous(breaks=c(1, 4, 8, 16), trans='sqrt') + 
  geom_point(aes(col=environment), size=2) + 
  geom_line(aes(col=environment), size = 0.5, alpha = 0.2) + 
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, col=environment), width=0.2) +
  facet_wrap(~ experiment, scales = 'free_y') +
  scale_color_grey() +
  xlab("Amount of computing units (count)") +
  ylab("Execution time (s)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", legend.spacing = unit(x=c(0,0,0,0),units="mm")) +
  default_theme()

plot(g)
