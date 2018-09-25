library(tidyverse)
library(RColorBrewer)

files <- list('./results/nas/results.csv',
              './results/ondes3d/results.csv',
              './results/ondes3d-ligurian/results.csv')

# Calculate makespan for all experiments
df.makespan <- do.call("bind_rows", lapply(files, function(file) {
  read_csv(file, col_types=cols(
    name = col_integer(),
    environment = col_character(),
    parallelism = col_integer(),
    time = col_integer()
  )) %>%
    mutate(time = time/1000) %>%
    group_by(environment, parallelism) %>%
    summarize(
      samples = n(),
      average = mean(time),
      stdDeviation = sd(time),
      stdError = 3*stdDeviation/sqrt(samples)
    ) %>%
    mutate(Origin = file) %>%
    separate(Origin, into=c("X0", "X1", "TYPE", "X2"), sep="/", remove=FALSE) %>% select(-X0, -X1, -X2) %>%
    mutate(Application = ifelse(grepl("ondes3d", TYPE), "Ondes3D", "NAS-EP")) %>%
    mutate(Input = case_when(TYPE == "nas" ~ "Class B",
                             grepl("ligurian", TYPE) ~ "Ligurian",
                             TRUE ~ "Default")) %>% select(-TYPE, -Origin) %>%
    mutate(Native = environment == "native") %>%
    ungroup()
})) %>%
  mutate(environment = factor(environment,
                              levels=c("native", "singularity", "docker"),
                              labels=c("Native", "Singularity", "Docker")))

# Create color palette
colors <- brewer.pal(3,"Set1")
myColors <- c(colors[1], colors[3], colors[2])
names(myColors) <- levels(df.sel$environment)
colScale <- scale_colour_manual(name = "environment",values = myColors)

# Plot RAW results...
df.makespan %>%
  mutate(Native = factor(Native, levels=c(TRUE, FALSE), labels=c("Native", "Container"))) %>%
  filter(Application == 'NAS-EP') %>%
  mutate(Application = paste(Application, Input, sep="/")) -> df.sel;

# Define breaks in X
breaks <- df.sel %>% pull(parallelism) %>% unique;

df.sel %>%
  ggplot(aes(x = parallelism, y = average, col=environment)) +
  ylim(0, NA) +
  geom_point(size=1) +
  geom_line(alpha = 0.2) + 
  geom_errorbar(aes(ymin = average - stdError, ymax = average + stdError, col = environment), width = 1) +
  colScale +
  xlab('Number of MPI ranks (count)') + 
  ylab('Execution time (s)') +
  theme(legend.position = 'top', legend.spacing = unit(x = c(0, 0, 0, 0), units = 'mm')) +
  theme_bw(base_size = 13) +
  theme (plot.margin = unit(c(0,0,0,0), "cm"),
         legend.spacing = unit(1, "mm"),
         legend.position = "top",
         legend.justification = "left",
         legend.box.spacing = unit(0, "pt"),
         legend.box.margin = margin(0,0,0,0),
         legend.title = element_blank()) +
  facet_grid(Application~Native, scales="free_y")

# Plot overhead of default cases
df.makespan.overhead <- df.makespan %>%
  mutate(Case = paste(Application, Input, sep="/")) %>%
  select(environment, parallelism, average, Case) %>%
  spread(environment, average) %>%
  mutate(Docker.Overhead = round((Docker - Native)/Native * 100, 2),
         Singularity.Overhead = round((Singularity - Native)/Native * 100, 2)) %>%
  arrange(Case, parallelism) %>%
  select(parallelism, Case, Docker.Overhead, Singularity.Overhead) %>%
  gather(environment, Overhead, -parallelism, -Case)

df.sel <- df.makespan.overhead %>%
  filter(Case == 'Ondes3D/Default') %>%
  mutate(environment = factor(gsub(".Overhead", "", environment),
                              levels=c("Singularity", "Docker")))
# Breaks in X
breaks <- df.sel %>% pull(parallelism) %>% unique;

df.sel %>%
  ggplot(aes(x=parallelism,
             y=Overhead,
             color=environment)) +
  scale_x_continuous(breaks = breaks, trans="sqrt") +
  colScale +
  geom_point(size=1) +
  geom_line(alpha = 0.3) +
  xlab('Number of MPI ranks (count)') + 
  ylab('Overhead of Execution \nTime against Native (%)') +
  facet_grid(Case~.) +
  ylim(0,60) +
  theme_bw(base_size = 13) +
  theme (plot.margin = unit(c(0,0,0,0), "cm"),
         legend.spacing = unit(x = c(0, 0, 0, 0), units = 'mm'),
         legend.position = "top",
         legend.justification = "left",
         legend.box.spacing = unit(0, "pt"),
         legend.box.margin = margin(0,0,0,0),
         legend.title = element_blank())

# Plot overhead of Ligurian
df.sel <- df.makespan.overhead %>%
  filter(Case == "Ondes3D/Ligurian") %>%
  mutate(environment = factor(gsub(".Overhead", "", environment),
                              levels=c("Singularity", "Docker"))) %>%
  filter(environment == "Singularity")

# Breaks in X
breaks <- df.sel %>% pull(parallelism) %>% unique;

df.sel %>%
  ggplot(aes(x=parallelism,
             y=Overhead,
             color=environment)) +
  scale_x_continuous(breaks = breaks, trans="sqrt") +
  colScale +
  geom_point(size=1) +
  geom_line(alpha = 0.3) +
  xlab('Number of MPI ranks (count)') + 
  ylab('Overhead of Execution \nTime against Native (%)') +
  facet_grid(Case~.) +
  ylim(0,10) +
  theme_bw(base_size = 13) +
  theme (plot.margin = unit(c(0,0,0,0), "cm"),
         legend.spacing = unit(x = c(0, 0, 0, 0), units = 'mm'),
         legend.position = "top",
         legend.justification = "left",
         legend.box.spacing = unit(0, "pt"),
         legend.box.margin = margin(0,0,0,0),
         legend.title = element_blank())
