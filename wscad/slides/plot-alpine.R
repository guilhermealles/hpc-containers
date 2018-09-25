library(tidyverse)
library(RColorBrewer)

# Read & mutate data
files <- list('./results/nas-alpine/results.csv',
              './results/ondes3d-alpine/results.csv')
df.alpine <- do.call("bind_rows", lapply(files, function(file) {
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
    mutate(Input = case_when(TYPE == "nas-alpine" ~ "Class B",
                             TYPE == "ondes3d-alpine" ~ "Default",
                             TRUE ~ "Default")) %>% select(-Origin) %>%
    mutate(Native = environment == "native") %>%
    ungroup()
})) %>%
  mutate(environment = factor(environment,
                              levels=c("native", "singularity", "docker"),
                              labels=c("Native", "Singularity", "Docker")))

# Plot
df.alpine %>%
  mutate(Native = factor(Native, levels=c(TRUE, FALSE), labels=c("Native", "Container"))) %>%
  mutate(Application = paste(Application, Input, sep="/")) %>%
  filter(Input %in% c("Class B", "Default")) -> df.sel;

# Create color palette
colors <- brewer.pal(3,"Set1")
myColors <- c(colors[1], colors[3], colors[2])
names(myColors) <- levels(df.sel$environment)
colScale <- scale_colour_manual(name = "environment",values = myColors)

# Breask in X
breaks <- df.sel %>% pull(parallelism) %>% unique;

# Plot raw results
df.sel %>%
  ggplot(aes(x = parallelism, y = average, col=environment)) +
  scale_x_continuous(breaks = breaks, trans="sqrt") +
  ylim(0, NA) +
  geom_point(size=1) +
  geom_line(alpha = 0.2) + 
  geom_errorbar(aes(ymin = average - stdError, ymax = average + stdError, col = environment), 
                width = .35) +
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
  facet_grid (Application~Native, scales="free_y")

# Plot overhead
df.alpine.overhead <- df.alpine %>%
  mutate(Case = paste(Application, Input, sep="/")) %>%
  select(environment, parallelism, average, Case) %>%
  spread(environment, average) %>%
  mutate(
    Docker.Gain = round((Docker - Native)/Native * 100, 2),
    Singularity.Gain = round((Singularity - Native)/Native * 100, 2)) %>%
  arrange(Case, parallelism) %>%
  select(parallelism, Case, Docker.Gain, Singularity.Gain) %>%
  gather(environment, Gain, -parallelism, -Case)

df.sel <- df.alpine.overhead %>%
  filter(parallelism <= 16) %>%
  mutate(environment = factor(gsub(".Gain", "", environment),
                              levels=c("Singularity", "Docker")))
# Breaks in X
breaks <- df.sel %>% pull(parallelism) %>% unique;

df.sel %>%
  ggplot(aes(x=parallelism,
             y=Gain,
             color=environment)) +
  scale_x_continuous(breaks = breaks, trans="sqrt") +
  colScale +
  geom_point(size=1) +
  geom_line(alpha = 0.3) +
  xlab('Number of MPI ranks (count)') + 
  ylab('Performance against Native (%)') +
  facet_grid(Case~.) +
  coord_cartesian(ylim=c(-20,20)) +
  #    ylim(-30,30) +
  theme_bw(base_size = 13) +
  theme (plot.margin = unit(c(0,0,0,0), "cm"),
         legend.spacing = unit(x = c(0, 0, 0, 0), units = 'mm'),
         legend.position = "top",
         legend.justification = "left",
         legend.box.spacing = unit(0, "pt"),
         legend.box.margin = margin(0,0,0,0),
         legend.title = element_blank())

