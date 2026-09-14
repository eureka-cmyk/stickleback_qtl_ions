library(tidyverse)

# input
pheno_file <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data.csv"
geno_file  <- "/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic/Sex_coverage/Sex_calling_summary.csv"

# output
out_file <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/sex_comparison.txt"

# phenotype sex
pheno <- read.csv(pheno_file) %>%
  select(fish_id, sex) %>%
  rename(ID = fish_id, pheno_sex = sex)

# genomic sex
geno <- read.csv(geno_file) %>%
  select(individual, sex) %>%
  rename(ID = individual, geno_sex = sex)

# combine
compare <- pheno %>%
  left_join(geno, by = "ID") %>%
  mutate(
    match = ifelse(pheno_sex == geno_sex, "YES", "NO")
  )

# save
write.table(compare,
            file = out_file,
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")


head(compare, 20)
