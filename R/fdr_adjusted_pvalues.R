library(tidyverse)
library(fdrtool)
par(mar=c(1,1,1,1))

Weights = as_tibble(atomX_hypothesis_3$W)

p_table = tibble(otu = unique(Weights$OTU))
p_table = p_table %>%
  mutate(weights = map(otu, ~filter(Weights, OTU == .))) %>%
  mutate(p_values = map(weights, ~t.test(W~IS_CONTROL, .)$p.value)) %>%
  mutate(p_values = as.numeric(p_values))

fdr_table = fdrtool(p_table$p_values, statistic = "pvalue")

p_table = p_table %>%
  mutate(adj_p_values = fdr_table$qval) %>%
  mutate(incorporator = adj_p_values < 0.05) %>%
  rename(OTU = otu)

incorporator_table = p_table %>% inner_join(df_atomX_boot_hypothesis_3 %>% select(-Incorporator))

ggplot(incorporator_table, aes(factor(OTU, levels = feature_factors,ordered = T), A, ymin=A_CI_low, ymax=A_CI_high, color=incorporator)) +
  geom_pointrange(size=0.25) +
  geom_linerange() +
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5) +
  labs(x='MAGs', y='Atom fraction excess') +
  theme_bw() +
  coord_flip() +
  ggtitle("Isotope incorporating MAGs") +
  ylim(0,0.6)

