# |----------------------------------------------------------------------------------|
# | Project: Human curcumin                                       |
# | Script: RNA-seq data analysis and visualization  |
# | Author: Davit Sargsyan                                                           |
# | Created: 02/25/2019                                                              |
# | Modified:                  |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_human_cur_v1.R.R")

require(data.table)
dt1 <- fread("data/human.curcumin.featurecounts.results_2.csv")
dt1

dt2 <- fread("data/human_cur_rna_sample_legend.csv")
dt2
