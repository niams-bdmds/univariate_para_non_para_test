library(dplyr) 
library(broom)
setwd("/path/to/work/dir")
data <- read.csv("countData.csv", header = TRUE, row.names = 1)

my_df <- data.frame(t(data)) %>% tibble::rownames_to_column(var = "gene-name")

sample_info <- read.table("sample_info.txt", header = TRUE)

data_new <- cbind(data, sample_info)

# Implement t test
t_test <- apply(data_new[1:(length(data_new)-5)], 2, function(i)tidy(t.test
        (i ~ data_new$Secondary_end_point))$p.value)
result_df <- bind_rows(t_test)
result_df <- t(result_df)
write.table(result_df, file = "Result_t_test_secondary_endpoint")

# Implement Mann-Whitney U test or Wilcox Rank Sum Test
wilcox_rank_test <- apply(data_new[1:(length(data_new)-5)], 2, function(i)tidy(wilcox.test
        (i ~ data_new$Primary_end_point))$p.value)
WRST_result <- bind_rows(wilcox_rank_test)
write.table(WRST_result, file = "Result_WRST_primary_end_point.txt")

# Two sample Kolmogorov-Smirnov test
KST <- apply(data_new[1:(length(data_new)-5)], 2, function(i)tidy(ks.test
        (i ~ data_new$Primary_end_point))$p.value)
KST_result <- bind_rows(KST)
write.table(KST_result, file = "Result_WRST_primary_end_point.txt")

rm(list=ls())

