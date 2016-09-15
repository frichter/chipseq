# Felix Richter
# 8/7/2016
# plot Chip-Seq enrichment

library("dplyr")
library("tidyr")
library("ggplot2")

setwd("D:/Dropbox/PhD/chipseq/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/chipseq")
setwd("/Users/felixrichter/Dropbox/PhD/chipseq")
options(stringsAsFactors=FALSE)


chip_enrich = read.table("pct_overlap.txt", sep = ",") %>% t() %>% as.data.frame
colnames(chip_enrich) = c("filename", "pct.overlap")
chip_enrich = chip_enrich %>%
  mutate(filename = gsub("input/", "", filename)) %>%
  mutate(filename = gsub(".merged", "", filename)) %>%
  separate(filename, into = c("new_chip", "db_anno"), sep = "___") %>%
  mutate(pct.overlap = as.numeric(pct.overlap) %>% signif(digits = 2))

p = chip_enrich %>%
  filter(!grepl("peaks", new_chip)) %>%
  mutate(db_anno = gsub("GSM.*UCSD.", "", db_anno)) %>% 
  mutate(db_anno = gsub(".*iPS-20b_H3K4Me3", "iPS-20b_H3K4Me3", db_anno)) %>% 
  mutate(db_anno = gsub(".*iPS-20b.H3K4Me3", "iPS-20b_H3K4Me3", db_anno)) %>% 
  mutate(db_anno = gsub("GSM.*UBC.", "", db_anno)) %>% 
  mutate(db_anno = gsub("GSM.*BI.", "", db_anno)) %>% 
  mutate(db_anno = gsub("regions_prom_E095", "LV_active_promoters", db_anno)) %>% 
  mutate(db_anno = gsub("regions_enh_E095", "LV_active_enhancers", db_anno)) %>% 
  mutate(db_anno = gsub("E095_15_coreMarks_1", "LV_TSS_ChromHMM", db_anno)) %>% 
  ggplot(., aes(x = new_chip, y = db_anno, fill = pct.overlap)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_bw() +
  geom_text(aes(label = pct.overlap), color="black") #+  #, size=1
  # theme(axis.text.x = element_text(angle = 30))
p
ggsave(filename = "pct_overlap.png", p, height = 5, width = 10, units = "in")
