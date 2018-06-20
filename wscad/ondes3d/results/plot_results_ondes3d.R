library(tidyverse);

docker_results <- read_csv('./doe-results.csv');
sing_native_results <- read_csv('./doe-results-sing-native.csv');

results <- bind_rows(sing_native_results, docker_results)
results <- results %>%
  mutate(time = time/1000) %>%
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
  geom_line(aes(col=environment), size = 0.5, alpha=0.2) + 
  geom_point(aes(col=environment), size=2) + 
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, col=environment), width=0.2) +
  scale_color_grey() +
  ylim(0, NA) +
  scale_x_continuous(breaks=c(1, 4, 8, 16), trans='sqrt') + 
  xlab("Amount of computing units (count)") +
  ylab("Execution time(s)") +
  theme_bw(base_size=12) +
  theme(legend.position = "top", legend.spacing = unit(x=c(0,0,0,0),units="mm")) +
  default_theme();

g2 <- g +
  xlab(NULL) +
  ylab(NULL) +
  coord_cartesian(ylim = c(3,10), xlim = c(4, 16)) +
  theme(legend.position = "none");

p2 <- ggplotGrob(g2);
g_pdf <- g + annotation_custom(grob = p2, xmin = 2, xmax = 4, ymin = 10, ymax = 32);

pdf(width=3.5, height=5, file='r-test-2.pdf');
plot(g_pdf);
dev.off()

