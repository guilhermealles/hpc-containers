library(tidyverse);

results <- read_csv('./doe-results.csv');
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
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, col=environment), width=20) +
  scale_color_grey() +
  scale_x_continuous(breaks=seq(64,256,64)) +
  ylim(0, NA) +
  xlab("Amount of computing units (count)") +
  ylab("Execution time (s)") +
  theme_bw(base_size=12) +
  theme(legend.position = "top", legend.spacing = unit(x=c(0,0,0,0),units="mm")) +
  default_theme();


pdf(width=3.5, height=5, file='r-test-2.pdf');
plot(g);
dev.off()

