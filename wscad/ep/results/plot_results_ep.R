library(tidyverse);

docker_results <- read_csv('./ep-results.csv');
sing_native_results <- read_csv('./ep-sing-native-results.csv');

results <- bind_rows(sing_native_results, docker_results);
results <- results %>%
  mutate(time=time/1000) %>%
  group_by(environment, parallelism) %>%
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
  ylim(0, NA) +
  geom_point(aes(col=environment), size=2) + 
  geom_line(aes(col=environment), size = 0.5, alpha = 0.2) + 
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, col=environment), width=0.2) +
  scale_color_grey() +
  xlab("Amount of computing units (count)") +
  ylab("Execution time (s)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", legend.spacing = unit(x=c(0,0,0,0),units="mm")) +
  default_theme();

g2 <- g +
  xlab(NULL) +
  ylab(NULL) +
  coord_cartesian(ylim = c(0,15), xlim = c(12, 20)) +
  theme(legend.position = "none");

p2 <- ggplotGrob(g2);
g3 <- g + annotation_custom(grob = p2, xmin = 2, xmax = 4, ymin = 30, ymax = 80);


plot(g3);

