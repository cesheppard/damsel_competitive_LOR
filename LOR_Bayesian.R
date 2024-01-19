library(tidyverse)
library(dplyr)
library(gridExtra)
library(boot)
library(cowplot)
library(sf)
library(ggpattern)
library(brms)
library(car)
library(jtools)
library(tidybayes)
library(emmeans)
library(posterior)
library(loo)
library(bayestestR)

### SET WORKING DIRECTORY ###
setwd("")

### READ IN DATA ###
damsel_behaviour <- read.csv("damsel_behaviour.csv", header = T, stringsAsFactors = F) # collated behavioural outputs from BORIS v.8.6.2 (Friard & Gamba, 2016) 
damsel_data <- read.csv("damsel_data.csv", header = T, stringsAsFactors = F) # individual damselfish data (species | transect | number of neighbours in buffers | territory size)
all_count_data <- read.csv("all_count_data.csv", header = T, stringsAsFactors = F) # buffer data (herbivore and non-herbivore counts | area | non-sand proportion) 

damsel_data <- damsel_data %>% 
  filter(id %in% damsel_behaviour$id) %>%
  arrange(desc(id)) %>%
  mutate(totalneigh1 = TSneigh1 + Dneigh1,
         totalneigh.5 = TSneigh.5 + Dneigh.5,
         totaldam1 = TSneigh1 + Dneigh1 + 1,
         totaldam.5 = TSneigh.5 + Dneigh.5 + 1) # total number of damsels within buffers

### CALCULATE AGGRESSION METRICS ###
damsel_aggression <- damsel_behaviour %>% 
  group_by(id) %>% 
  summarise(total_intruders = n())

for (i in 1:length(damsel_behaviour$id)) {
  damsel_behaviour$dietary[i] <- ifelse(damsel_behaviour$family[i] %in% c("Chub", "Damselfish", "Parrotfish","Pufferfish", "Surgeonfish"), "Herbivore", "None")
} # label herbivores in behavioural data

damsel_aggression <- left_join(damsel_behaviour %>%
                                 group_by(id) %>%
                                 summarise(total_intruders = n()),
                               damsel_behaviour %>%
                                 group_by(id) %>%
                                 filter(behaviour == "aggression") %>%
                                 summarise(total_chases = n())) %>%
  left_join(damsel_behaviour %>%
              group_by(id) %>%
              filter(type == "Heterospecific" | type == "None") %>%
              summarise(hetero_intruders = n())) %>%
  left_join(damsel_behaviour %>%
              group_by(id) %>%
              filter(type == "Heterospecific" & behaviour == "aggression" | type == "None" & behaviour == "aggression") %>%
              summarise(hetero_chases = n())) %>%
  left_join(damsel_behaviour %>%
              group_by(id) %>%
              filter(type == "Heterospecific" & dietary == "Herbivore" |
                       type == "None" & dietary == "Herbivore") %>%
              summarise(herb_intruders = n())) %>%
  left_join(damsel_behaviour %>%
              group_by(id) %>%
              filter(type == "Heterospecific" & dietary == "Herbivore" & behaviour == "aggression" |
                       type == "None" & dietary == "Herbivore" & behaviour == "aggression" ) %>%
              summarise(herb_chases = n())) %>%
  left_join(damsel_behaviour %>%
              group_by(id) %>%
              filter(type == "Heterospecific" & dietary == "None" |
                       type == "None" & dietary == "None") %>%
              summarise(non_herb_intruders = n())) %>%
  left_join(damsel_behaviour %>%
              group_by(id) %>%
              filter(type == "Heterospecific" & dietary == "None" & behaviour == "aggression" |
                       type == "None" & dietary == "None" & behaviour == "aggression") %>%
              summarise(non_herb_chases = n()))

damsel_aggression <- left_join(damsel_aggression, damsel_data)
damsel_aggression <- damsel_aggression %>% 
  mutate_at(c(2:9), ~replace_na(., 0))

damsel_aggression <- damsel_aggression %>%
  mutate(freq = total_chases/total_intruders,
         hetero.freq = hetero_chases/hetero_intruders,
         herb.freq = herb_chases/herb_intruders,
         non.herb.freq = non_herb_chases/non_herb_intruders) # aggression metrics based on proportion of intruders chased

### BUFFER COUNTS AND AREAS ###
final <- damsel_aggression %>%
  left_join(all_count_data) %>%
  mutate_all(~replace_na(., 0))
final$transect <- as.factor(final$transect)
# binomial yes|no of whether there is a sand patch within the buffer also explored as many buffers have no sand

### BAYESIAN MODELLING ###
# total count 1 m buffer
priors <- c(set_prior("normal(0, 10)", class = "b")) # weakly informative

total1.model <- brm(totalcount1 ~ hetero.freq + totaldam1 + (1|transect) + offset(log(area1)),
                    prior = priors, 
                    family = negbinomial(),
                    iter = 5000, warmup = 1000, chains = 4, cores = 4, 
                    control=list(adapt_delta=0.95),
                    data = final)
summary(total1.model)

# check plots
pp_check(total1.model, ndraws = 10) # looks okay
plot(total1.model, ask = FALSE)
bayes_R2(total1.model)

# check for influential data points
total1.loo <- loo(total1.model)
pareto_k_ids(total1.loo, threshold = .7)

# model coefficients
p_direction(total1.model)

aggression_coef <- as.data.frame(fixef(total1.model, summary = FALSE) [,2])
names(aggression_coef) <- c("coef")
aggression_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") # visualize coefficients

damsel_coef <- as.data.frame(fixef(total1.model, summary = FALSE) [,3])
names(damsel_coef) <- c("coef")
damsel_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed")

# make plot
plot_theme <-
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        strip.background = element_blank(), 
        strip.text = element_text(size = 10))

facet.labs <- c("Transect 1", "Transect 2")
names(facet.labs) <- c("1", "2")

pred_total1 <-
  total1.model %>% 
  epred_draws(newdata = expand_grid(hetero.freq = seq(from = 0, to = 1, by = 0.01),
                                    transect =  c("1", "2"),
                                    totaldam1 = mean(final$totaldam1),
                                    area1 = mean(final$area1),
                                    reps = 100), 
              re_formula=NULL)

plot1 <- ggplot(aes(x = hetero.freq, y = totalcount1), data = final) +
  geom_point() +
  stat_lineribbon(aes(y = .epred), data = pred_total1, .width = c(.95, .80, .50), alpha = .5) +
  scale_fill_brewer(palette = "Blues") +
  #scale_x_continuous(name = "Aggression (focal Stegastes towards non-Stegastes intruders)") +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "1 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2, labeller = labeller(transect = facet.labs)) +
  plot_theme
plot1

# total count .5 m buffer
total5.model <- brm(totalcount.5 ~ hetero.freq + totaldam.5 + (1|transect) + offset(log(area.5)),
                    prior = priors, 
                    family = negbinomial(),
                    iter = 5000, warmup = 1000, chains = 4, cores = 4, 
                    control=list(adapt_delta=0.95),
                    data = final)
summary(total5.model)

# check plots
pp_check(total5.model, ndraws = 10) # looks okay
plot(total5.model, ask = FALSE)
bayes_R2(total5.model)

# check for influential data points
total5.loo <- loo(total5.model)
pareto_k_ids(total5.loo, threshold = .7)

# model coefficients
p_direction(total5.model)

aggression_coef <- as.data.frame(fixef(total5.model, summary = FALSE) [,2])
names(aggression_coef) <- c("coef")
aggression_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") # visualize coefficients

damsel_coef <- as.data.frame(fixef(total5.model, summary = FALSE) [,3])
names(damsel_coef) <- c("coef")
damsel_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed")

# aggression plot
pred_total5 <-
  total5.model %>% 
  epred_draws(newdata = expand_grid(hetero.freq = seq(from = 0, to = 1, by = 0.01),
                                    transect =  c("1", "2"),
                                    totaldam.5 = mean(final$totaldam.5),
                                    area.5 = mean(final$area.5),
                                    reps = 100), 
              re_formula=NULL)

plot2 <- ggplot(aes(x = hetero.freq, y = totalcount.5), data = final) +
  geom_point() +
  stat_lineribbon(aes(y = .epred), data = pred_total5, .width = c(.95, .80, .50), alpha = .5) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(name = "Aggression (focal Stegastes towards non-Stegastes intruders)") +
  scale_y_continuous(name = "0.5 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2) +
  plot_theme +
  theme(strip.text = element_blank())
plot2

# herb count 1 m buffer
herb1.model <- brm(count1 ~ herb.freq + totaldam1 + (1|transect) + offset(log(area1)),
                    prior = priors, 
                    family = negbinomial(),
                    iter = 5000, warmup = 1000, chains = 4, cores = 4, 
                    control=list(adapt_delta=0.95),
                    data = final)
summary(herb1.model)

# check plots
pp_check(herb1.model, ndraws = 10) # looks okay
plot(herb1.model, ask = FALSE)
bayes_R2(herb1.model)

# check for influential data points
herb1.loo <- loo(herb1.model)
pareto_k_ids(herb1.loo, threshold = .7)

# model coefficients
p_direction(herb1.model)

aggression_coef <- as.data.frame(fixef(herb1.model, summary = FALSE) [,2])
names(aggression_coef) <- c("coef")
aggression_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") # visualize coefficients

damsel_coef <- as.data.frame(fixef(herb1.model, summary = FALSE) [,3])
names(damsel_coef) <- c("coef")
damsel_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed")

# aggression plot
pred_herb1 <-
  herb1.model %>% 
  epred_draws(newdata = expand_grid(herb.freq = seq(from = 0, to = 1, by = 0.01),
                                    transect =  c("1", "2"),
                                    totaldam1 = mean(final$totaldam1),
                                    area1 = mean(final$area1),
                                    reps = 100), 
              re_formula=NULL)

plot3 <- ggplot(aes(x = herb.freq, y = count1), data = final) +
  geom_point() +
  stat_lineribbon(aes(y = .epred), data = pred_herb1, .width = c(.95, .80, .50), alpha = .5) +
  scale_fill_brewer(palette = "Blues") +
  #scale_x_continuous(name = "Aggression (focal Stegastes towards herbivorous intruders)") +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "1 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2, labeller = labeller(transect = facet.labs)) +
  plot_theme
plot3

# non herb count 1 m buffer
nonherb1.model <- brm(nhcount1 ~ non.herb.freq + totaldam1 + (1|transect) + offset(log(area1)),
                   prior = priors, 
                   family = negbinomial(),
                   iter = 5000, warmup = 1000, chains = 4, cores = 4, 
                   control=list(adapt_delta=0.95),
                   data = final)
summary(nonherb1.model)

# check plots
pp_check(nonherb1.model, ndraws = 10) # seems okay
plot(nonherb1.model, ask = FALSE)
bayes_R2(nonherb1.model)

# check for influential data points
nonherb1.loo <- loo(nonherb1.model)
pareto_k_ids(nonherb1.loo, threshold = .7)

# model coefficients
p_direction(nonherb1.model)

aggression_coef <- as.data.frame(fixef(nonherb1.model, summary = FALSE) [,2])
names(aggression_coef) <- c("coef")
aggression_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") # visualize coefficients

damsel_coef <- as.data.frame(fixef(nonherb1.model, summary = FALSE) [,3])
names(damsel_coef) <- c("coef")
damsel_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed")

# aggression plot
pred_nonherb1 <-
  nonherb1.model %>% 
  epred_draws(newdata = expand_grid(non.herb.freq = seq(from = 0, to = 1, by = 0.01),
                                    transect =  c("1", "2"),
                                    totaldam1 = mean(final$totaldam1),
                                    area1 = mean(final$area1),
                                    reps = 100), 
              re_formula=NULL)

plot4 <- ggplot(aes(x = non.herb.freq, y = nhcount1), data = final) +
  geom_point() +
  stat_lineribbon(aes(y = .epred), data = pred_nonherb1, .width = c(.95, .80, .50), alpha = .5) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(name = "Aggression (focal Stegastes towards non herbivorous intruders)") +
  scale_y_continuous(name = "1 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2, labeller = labeller(transect = facet.labs)) +
  plot_theme
plot4

# herb count .5 m buffer
herb5.model <- brm(count.5 ~ herb.freq + totaldam.5 + (1|transect) + offset(log(area.5)),
                   prior = priors, 
                   family = negbinomial(),
                   iter = 5000, warmup = 1000, chains = 4, cores = 4, 
                   control=list(adapt_delta=0.95),
                   data = final)
summary(herb5.model)

# check plots
pp_check(herb5.model, ndraws = 10) # looks okay
plot(herb5.model, ask = FALSE)
bayes_R2(herb5.model)

# check for influential data points
herb5.loo <- loo(herb5.model)
pareto_k_ids(herb5.loo, threshold = .7)

# model coefficients
p_direction(herb5.model)

aggression_coef <- as.data.frame(fixef(herb5.model, summary = FALSE) [,2])
names(aggression_coef) <- c("coef")
aggression_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") # visualize coefficients

damsel_coef <- as.data.frame(fixef(herb5.model, summary = FALSE) [,3])
names(damsel_coef) <- c("coef")
damsel_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed")

# aggression plot
pred_herb5 <-
  herb5.model %>% 
  epred_draws(newdata = expand_grid(herb.freq = seq(from = 0, to = 1, by = 0.01),
                                    transect =  c("1", "2"),
                                    totaldam.5 = mean(final$totaldam.5),
                                    area.5 = mean(final$area.5),
                                    reps = 100), 
              re_formula=NULL)

plot5 <- ggplot(aes(x = herb.freq, y = count.5), data = final) +
  geom_point() +
  stat_lineribbon(aes(y = .epred), data = pred_herb5, .width = c(.95, .80, .50), alpha = .5) +
  scale_fill_brewer(palette = "Blues") +
  #scale_x_continuous(name = "Aggression (focal Stegastes towards herbivorous intruders)") +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "0.5 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2) +
  plot_theme +
  theme(strip.text = element_blank())
plot5

# non herb count .5 m buffer
nonherb5.model <- brm(nhcount.5 ~ non.herb.freq + totaldam.5 + (1|transect) + offset(log(area.5)),
                      prior = priors, 
                      family = negbinomial(),
                      iter = 5000, warmup = 1000, chains = 4, cores = 4, 
                      control=list(adapt_delta=0.95),
                      data = final)
summary(nonherb5.model)

# check plots
pp_check(nonherb5.model, ndraws = 10) # looks okay
plot(nonherb5.model, ask = FALSE)
bayes_R2(nonherb5.model)

# check for influential data points
nonherb5.loo <- loo(nonherb5.model)
pareto_k_ids(nonherb5.loo, threshold = .7)

# model coefficients
p_direction(nonherb5.model)

aggression_coef <- as.data.frame(fixef(nonherb5.model, summary = FALSE) [,2])
names(aggression_coef) <- c("coef")
aggression_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed") # visualize coefficients

damsel_coef <- as.data.frame(fixef(nonherb5.model, summary = FALSE) [,3])
names(damsel_coef) <- c("coef")
damsel_coef %>% 
  ggplot(aes(x = coef)) +
  stat_slab(alpha = .5) +
  geom_vline(xintercept = 0, linetype = "dashed")

# aggression plot
pred_nonherb5 <-
  nonherb5.model %>% 
  epred_draws(newdata = expand_grid(non.herb.freq = seq(from = 0, to = 1, by = 0.01),
                                    transect =  c("1", "2"),
                                    totaldam.5 = mean(final$totaldam.5),
                                    area.5 = mean(final$area.5),
                                    reps = 100), 
              re_formula=NULL)

plot6 <- ggplot(aes(x = non.herb.freq, y = nhcount.5), data = final) +
  geom_point() +
  stat_lineribbon(aes(y = .epred), data = pred_nonherb5, .width = c(.95, .80, .50), alpha = .5) +
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(name = "Aggression (focal Stegastes towards non herbivorous intruders)") +
  scale_y_continuous(name = "0.5 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2) +
  plot_theme +
  theme(strip.text = element_blank())
plot6

# manuscript figures
# figure 2
plot_grid(plot1, plot2, ncol=1, align = "vh", axis = "rlbt", scale = .95) +
  draw_label("Non-Stegastes fish in surrounding area", x = 0, y = .5, vjust = 2, angle = 90, size = 10)

# figure S2
plot_grid(plot3, plot5, ncol=1, align = "vh", axis = "rlbt", scale = .95) +
  draw_label("Non-Stegastes herbivorous fish in surrounding area", x = 0, y = .5, vjust = 2, angle = 90, size = 10)

# figure S3
plot_grid(plot4, plot6, ncol=1, align = "vh", axis = "rlbt", scale = .95) +
  draw_label("Non-Stegastes non-herbivorous fish in surrounding area", x = 0, y = .5, vjust = 2, angle = 90, size = 10)


# damsel abundance plots
# total count 1 m buffer
pred_total1dams <-
  total1.model %>% 
  epred_draws(newdata = expand_grid(hetero.freq = mean(final$hetero.freq),
                                    transect =  c("1", "2"),
                                    totaldam1 = seq(from = 2, to = 6, by = 1),
                                    area1 = mean(final$area1),
                                    reps = 100), 
              re_formula=NULL) %>%
  filter(!(transect == 2 & totaldam1 == 6))


box1 <- ggplot() +
  geom_point(aes(x = totaldam1 - .1, y = totalcount1), data = final, colour = "#6baed6") +
  stat_pointinterval(aes(x = totaldam1 + .1, y = .epred), data = pred_total1dams,
                     point_interval = median_hdi, .width = c(.95, .80), fatten_point = 2, 
                     interval_alpha = 1, position = position_dodge(width = .2)) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "1 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2, labeller = labeller(transect = facet.labs)) +
  plot_theme
box1

# total count .5 m buffer
pred_total5dams <-
  total5.model %>% 
  epred_draws(newdata = expand_grid(hetero.freq = mean(final$hetero.freq),
                                    transect =  c("1", "2"),
                                    totaldam.5 = seq(from = 2, to = 4, by = 1),
                                    area.5 = mean(final$area.5),
                                    reps = 100), 
              re_formula=NULL)

box2 <- ggplot() +
  geom_point(aes(x = totaldam.5 - .1, y = totalcount.5), data = final, colour = "#6baed6") +
  stat_pointinterval(aes(x = totaldam.5 + .1, y = .epred), data = pred_total5dams,
                     point_interval = median_hdi, .width = c(.95, .80), fatten_point = 2, 
                     interval_alpha = 1, position = position_dodge(width = .2)) +
  scale_x_continuous(name = "Stegastes spp. abundance") +
  scale_y_continuous(name = "0.5 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2) +
  plot_theme +
  theme(strip.text = element_blank())
box2

# herb count 1 m buffer
pred_herb1dams <-
  herb1.model %>% 
  epred_draws(newdata = expand_grid(herb.freq = mean(final$herb.freq),
                                    transect =  c("1", "2"),
                                    totaldam1 = seq(from = 3, to = 6, by = 1),
                                    area1 = mean(final$area1),
                                    reps = 100), 
              re_formula=NULL) %>%
  filter(!(transect == 2 & totaldam1 == 6))


box3 <- ggplot() +
  geom_point(aes(x = totaldam1 - .1, y = count1), data = final, colour = "#6baed6") +
  stat_pointinterval(aes(x = totaldam1 + .1, y = .epred), data = pred_herb1dams,
                     point_interval = median_hdi, .width = c(.95, .80), fatten_point = 2, 
                     interval_alpha = 1, position = position_dodge(width = .2)) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "1 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2, labeller = labeller(transect = facet.labs)) +
  plot_theme
box3

# non herb count 1 m buffer
pred_nonherb1dams <-
  nonherb1.model %>% 
  epred_draws(newdata = expand_grid(non.herb.freq = mean(final$non.herb.freq),
                                    transect =  c("1", "2"),
                                    totaldam1 = seq(from = 3, to = 6, by = 1),
                                    area1 = mean(final$area1),
                                    reps = 100), 
              re_formula=NULL) %>%
  filter(!(transect == 2 & totaldam1 == 6))


box4 <- ggplot() +
  geom_point(aes(x = totaldam1 - .1, y = nhcount1), data = final, colour = "#6baed6") +
  stat_pointinterval(aes(x = totaldam1 + .1, y = .epred), data = pred_nonherb1dams,
                     point_interval = median_hdi, .width = c(.95, .80), fatten_point = 2, 
                     interval_alpha = 1, position = position_dodge(width = .2)) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "1 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2, labeller = labeller(transect = facet.labs)) +
  plot_theme
box4

# herb count .5 m buffer
pred_herb5dams <-
  herb5.model %>% 
  epred_draws(newdata = expand_grid(herb.freq = mean(final$herb.freq),
                                    transect =  c("1", "2"),
                                    totaldam.5 = seq(from = 2, to = 4, by = 1),
                                    area.5 = mean(final$area.5),
                                    reps = 100), 
              re_formula=NULL)

box5 <- ggplot() +
  geom_point(aes(x = totaldam.5 - .1, y = count.5), data = final, colour = "#6baed6") +
  stat_pointinterval(aes(x = totaldam.5 + .1, y = .epred), data = pred_herb5dams,
                     point_interval = median_hdi, .width = c(.95, .80), fatten_point = 2, 
                     interval_alpha = 1, position = position_dodge(width = .2)) +
  scale_x_continuous(name = "Stegastes spp. abundance") +
  scale_y_continuous(name = "0.5 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2) +
  plot_theme +
  theme(strip.text = element_blank())
box5

# non herb count .5 m buffer
pred_nonherb5dams <-
  nonherb5.model %>% 
  epred_draws(newdata = expand_grid(non.herb.freq = mean(final$non.herb.freq),
                                    transect =  c("1", "2"),
                                    totaldam.5 = seq(from = 2, to = 4, by = 1),
                                    area.5 = mean(final$area.5),
                                    reps = 100), 
              re_formula=NULL)

box6 <- ggplot() +
  geom_point(aes(x = totaldam.5 - .1, y = nhcount.5), data = final, colour = "#6baed6") +
  stat_pointinterval(aes(x = totaldam.5 + .1, y = .epred), data = pred_nonherb5dams,
                     point_interval = median_hdi, .width = c(.95, .80), fatten_point = 2, 
                     interval_alpha = 1, position = position_dodge(width = .2)) +
  scale_x_continuous(name = "Stegastes spp. abundance") +
  scale_y_continuous(name = "0.5 m buffer") +
  facet_wrap(~ transect, scales = "free", ncol = 2) +
  plot_theme +
  theme(strip.text = element_blank())
box6

# manuscript figures
# figure 2
plot_grid(box1, box2, ncol=1, align = "vh", axis = "rlbt", scale = .95) +
  draw_label("Non-Stegastes fish in surrounding area", x = 0, y = .5, vjust = 2, angle = 90, size = 10)

# figure S2
plot_grid(box4, box6, ncol=1, align = "vh", axis = "rlbt", scale = .95) +
  draw_label("Non-Stegastes non-herbivorous fish in surrounding area", x = 0, y = .5, vjust = 2, angle = 90, size = 10)

# figure S3
plot_grid(box3, box5, ncol=1, align = "vh", axis = "rlbt", scale = .95) +
  draw_label("Non-Stegastes herbivorous fish in surrounding area", x = 0, y = .5, vjust = 2, angle = 90, size = 10)
