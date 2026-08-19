library(tidyverse)

reads <- read_csv("data/read_amounts.csv", col_names = c("File","Reads"))

sum(reads$Reads)
controls <- c("KH_D024", "KH_D048", "KH_D049","KH_D050","KH_D084","KH_D085","KH_D086", "ZYMO_MOCK","PCR_BLANK", "F12_Blank06")

reads %>%
  separate_wider_delim(File, delim = "_", names = c("...1","...2","...3","x","y","KH","Stn","z","d", "e"), too_few = "align_end") %>%
  select(KH, Stn, Reads) %>% 
  mutate(ID = paste(KH, Stn, sep = "_")) %>%
  filter(!ID %in% controls) %>% 
  summarize(median(Reads),
            min(Reads),
            max(Reads))

abs_physeq
abs_physeq

sample_sums(floor_physeq) %>% sort()

sample_sums(euk_diversity_physeq) %>% sort()
