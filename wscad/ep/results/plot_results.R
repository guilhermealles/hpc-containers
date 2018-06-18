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
    stdError = 2*stdDeviation/sqrt(samples)
  );

g <- ggplot(results, aes(x = parallelism, y = average)) + 
  scale_x_continuous(breaks=c(1, 4, 8, 16), trans='sqrt') + 
  scale_y_continuous(breaks=seq(0, 90, 10), trans='sqrt') +
  geom_point(aes(col=environment), size=2) + 
  geom_line(aes(col=environment), size = 0.5) + 
  geom_errorbar(aes(ymin=average-stdError, ymax=average+stdError, col=environment), width=0.5) +
  ggtitle("NAS EP Benchmark", subtitle = "Category B") +
  xlab("MPI Processes") +
  ylab("Execution time(s)");

plot(g);
