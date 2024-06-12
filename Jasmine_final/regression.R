library(LaplacesDemon)
library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
library(randomcoloR)
library('lmtest')
setwd('~/Documents/qfm2')
devtools::load_all()

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()

setwd('~/Documents/R/Jasmine/')
indelphi_table <- readr::read_tsv('MARC1_indelphi_result.txt', col_names = T)
mut_rate_data <- read_csv('SupplementaryTable3 - TableS3-hgRNA Activities.csv')


raw_data = load_samples('~/Downloads/Jasmine Data (1)/MARC1-Pipeline-MARCC-master/5-pair_filtering/')
parent = load_parent('./PB3-founder_filteredpairs.txt')
raw_data = id_filter(raw_data, parent)
raw_data = annotate_parental_spacers(raw_data, parent)
raw_data = compute_mosaic_fraction(raw_data)
sample_metadata = as_tibble(stringr::str_match(unique(raw_data$data$sample), "(t\\d+)-(\\d+)-(.*?)_filteredpairs\\.txt"))
colnames(sample_metadata) = c('sample', 'mother', 'mouse', 'tissue')
sample_metadata$condition = map_chr(sample_metadata$mother, function(m) {
  if (grepl('t67', m)) {
    'MARC1'
  } else {
    'Cas9'
  }
})
raw_data = append_metadata(raw_data, sample_metadata)
colnames(mut_rate_data)[3] = 'ID'
raw_data = append_metadata(raw_data, mut_rate_data, cols = c('Class'), join_by = 'ID')
raw_data$data$tissue = substr(raw_data$data$tissue, 1, 2)
raw_data$data = raw_data$data[raw_data$data$tissue != 'MA' &
                              raw_data$data$tissue != 'pl',]

Cas9_IDs = unique(raw_data$data$ID[raw_data$data$condition == 'Cas9'])
MARC1_IDs = unique(raw_data$data$ID[raw_data$data$condition != 'Cas9'])
intersect(Cas9_IDs, MARC1_IDs)

#raw_data$data = raw_data$data %>% filter(Class != 'inactive')
colnames(indelphi_table)[c(2,4)] = c('sequence', 'probability')
#raw_data = filter_recurring_spacers(raw_data, indelphi_table, 0.001, filter_missing = FALSE)
nested_data = nest_raw(raw_data, c('ID', 'sequence'))
num_tissue_condition = map_dbl(nested_data$nested$data, function(data) {
  nrow(unique(select(ungroup(data), condition, tissue)))
})
nested_data$nested = nested_data$nested[num_tissue_condition == 6, ]
raw_data$data = unnest(nested_data$nested, cols = data)
nested_data = nest_raw(raw_data, c('ID', 'sample', 'tissue', 'condition', 'Class'))
nested_data = compute_mutation_information(nested_data)
# nested_data$nested$mutated_fraction = map_dbl(nested_data$nested$data, function(data) {
#   sum(data$mosaic_fraction)
# })

tb = nested_data$nested %>% group_by(condition, ID) %>%
  summarise(MF = mean(mutated_fraction),
            sd = sd(mutated_fraction))
plot_tb = tibble(Cas9 = tb$MF[tb$condition == 'Cas9'],
                 MARC1 = tb$MF[tb$condition != 'Cas9'],
                 sdCas9 = tb$sd[tb$condition == 'Cas9'],
                 sdMARC1 = tb$sd[tb$condition != 'Cas9'])
plot_tb = plot_tb %>% filter(Cas9 > 0)
plot_tb$ones = (1:nrow(plot_tb))/nrow(plot_tb)
plot_tb$Cas9 = jitter(plot_tb$Cas9, 20)
plot_tb$MARC1 = jitter(plot_tb$MARC1, 20)
dotplot = ggplot(data = plot_tb,
       aes(x = Cas9,
           y = MARC1)) +
  geom_point() +
  geom_errorbar(aes(y = MARC1,
                    xmin = Cas9 - sdCas9,
                    xmax = Cas9 + sdCas9)) +
  geom_errorbar(aes(x = Cas9,
                    ymin = MARC1 - sdMARC1,
                    ymax = MARC1 + sdMARC1)) +
  xlab('Cas9 Female x MARC1 Male') +
  ylab('MARC1 Female x Cas9 Male') +
  geom_line(aes(x = ones,
                y = ones))

png(filename = '~/Downloads/dotplot.png', width = 4, height = 4, units = "in", res = 500)
print(dotplot)
dev.off()


#bar plots
x = nested_data$nested %>% nest(n_data = -c(tissue, condition, Class)) %>%
  mutate(coef = map(n_data, function(tb) {
    y = glm(data = tb, formula = mutated_fraction ~ 1, family = "binomial") %>% 
      summary()
    tibble(mf = inverse_logit(y$coefficients[1,1]),
           mf_lo = inverse_logit(y$coefficients[1,1] - y$coefficients[1,2]),
           mf_hi = inverse_logit(y$coefficients[1,1] + y$coefficients[1,2]))
  })) %>% select(-n_data) %>%
  unnest(cols = coef)

barplot = ggplot(data = x) +
  geom_bar(aes(x = condition, y = mf, fill = tissue), position = position_dodge(),stat="identity") +
  geom_errorbar(aes(x = condition,
                    ymin = mf_lo,
                    ymax = mf_hi,
                    color = tissue),
                position = position_dodge(),
                stat="identity") +
  facet_wrap(~Class) +
  ylab('Mutated Fraction') +
  xlab('Mother Background') +
  theme(legend.position = 'right',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9))


inverse_logit(2.9687)
inverse_logit(2.9687 + 0.6993)
inverse_logit(2.9687 - 0.6993)


plot_tb = nested_data$nested %>% group_by(tissue, condition, Class) %>%
  summarise(MF = inverse_logit(mean(logit(mutated_fraction))),
            upper = mean(logit(mutated_fraction)) + se(logit(mutated_fraction)),
            lower = mean(logit(mutated_fraction)) - se(logit(mutated_fraction)))
ggplot(data = plot_tb) +
  geom_bar(aes(x = tissue, y = MF, fill = condition), stat="identity") +
  geom_errorbar(aes(x = tissue,
                    ymin = inverse_logit(lower),
                    ymax = inverse_logit(upper),
                    color = condition), stat="identity") +
  facet_wrap(~Class) +
  ylab('Mutated Fraction')

barplot = ggplot(data = plot_tb, aes(x = condition, y = MF, fill = tissue)) +
  scale_fill_manual(values = c('he' = 'red',
                               'ta' = 'blue',
                               'yo' = 'orange')) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(x = condition,
                    ymin = inverse_logit(lower),
                    ymax = inverse_logit(upper)),
                position = position_dodge()) +
  facet_wrap(~Class) + 
  theme(legend.position = 'right',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 9)) +
  ylab('Mutated Fraction') +
  xlab('Mother Background')

png(filename = '~/Downloads/color_tissue.png', width = 4, height = 4, units = "in", res = 500)
print(barplot)
dev.off()

#for all IDs combined

logit <- function(p) {
  log((p + 10^-10)/(1 - p + 10^-10))
}

inverse_logit <- function(x) {
  1/(1+exp(-x))
}

se <- function(x) {
  sd(x)/sqrt(length(x))
}

model = lm(formula = logit(mutated_fraction) ~ tissue + condition, data = nested_data$nested)
summary(model)
model_with_interaction_terms = lm(formula = logit(mutated_fraction) ~ tissue * condition,
                                  data = nested_data$nested)
summary(model_with_interaction_terms)
lrtest(model, model_with_interaction_terms)

#stratified by ID
nested_data_id = nest(nested_data$nested, data2 = -Class)
#model_1 = lm(formula = log(mutated_fraction + 0.0001) ~ tissue + condition, data = nested_data_id$data2[where(nested_data_id$ID == 'ATCG')])
nested_data_id$Class[[2]]

nested_data_id$model = map(nested_data_id$data2, function(data2) {
  lm(formula = logit(mutated_fraction) ~ tissue + condition, data = data2)
})
nested_data_id$model_with_interaction = map(nested_data_id$data2, function(data2) {
  lm(formula = logit(mutated_fraction) ~ tissue * condition, data = data2)
})




