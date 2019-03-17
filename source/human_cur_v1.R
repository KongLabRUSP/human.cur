# |----------------------------------------------------------------------------------|
# | Project: Human curcumin                                                          |
# | Script: RNA-seq data analysis and visualization                                  |
# | Author: Davit Sargsyan                                                           |
# | Created: 02/25/2019                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_human_cur_v1.R.R")

# Install packages----
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")
# BiocManager::install("org.Mm.eg.db", version = "3.8")
# BiocManager::install("BiocParallel", version = "3.8")

# Load packages----
require(data.table)
require(plotly)
# require(org.Mm.eg.db)
# require(DESeq2)
# require(BiocParallel)

# Load RNA counts----
dt1 <- fread("data/human.curcumin.featurecounts.results_2.csv")
dt1
dt1 <- data.table(gene = dt1$Geneid,
                  dt1[, Length:`S9-4h.sorted.dedup.bam`])
colnames(dt1)[-c(1:2)] <- gsub(pattern = '.sorted.dedup.bam',
                          replacement = "",
                          x = colnames(dt1)[-c(1:2)])
dt1
length(unique(dt1$gene))
# 26,364, same as the number of rows

# NOTE: there are 2 S7-0h samples! Removing `S7-0h-to`!
dt1 <- dt1[, -c("S7-0h-to")]

# Reorder columns----
dt1 <- dt1[, c(1:2, 5:14, 3:4)]
dt1

# Number of unique reads in each samle----
nreads <- data.table(Sample = colnames(dt1)[-1],
                     `Number of Total Hits` = colSums(dt1[, -1]))
nreads

# Load sample legend----
map <- fread("data/human_cur_rna_sample_legend.csv")
map <- map[, c(3:5, 2, 6)]
colnames(map)[1] <- "Sample"
map

# Merge----
nreads$Sample
map$Sample
nreads <- merge(nreads,
                map,
                by = "Sample")
nreads$Sample <- factor(nreads$Sample,
                        levels = unique(map$Sample))
nreads

# Plot number of unique reads per sample----
p1 <- ggplot(nreads,
             aes(x = Sample,
                 y = `Number of Total Hits`,
                 fill = Sex,
                 group = `Conc (ng/ÂµL)`)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# Plot correlation of total hits and RNA concentration----
p2 <- ggplot(nreads,
             aes(x = `Conc (ng/ÂµL)`,
                 y = `Number of Total Hits`,
                 fill = Sex,
                 group = Sample)) +
  geom_point(size = 3,
             alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p2)
# No correlation between the sequencing depth and the sample amount

# Relative abundane of genes per sample----
# NOTE: I used Transcripts Per Kilobase Million (TMP) normalization here, which is, for the j-th sample:
# 1. normilize for gene length: a[i, j] = 1,000*count[i, j]/gene[i, j] length(bp)
# 2. normalize for seq depth (i.e. total count): a(i, j)/sum(a[, j])
# 3. multiply by one million
# Follow links in the README file for a very good comparison of normalization techniques
tmp <- 1000*dt1[, 3:14]/dt1$Length
tmp
ra <- data.table(gene = dt1$gene,
                 apply(tmp,
                       2,
                       function(a) {
                         10^6*(a/sum(a))
                       }))
ra
colSums(ra[, -1])
# All 1e+06; CORRECT

summary(ra[, -1])
hist(log(as.matrix(ra[, -1] + 1)), 100)
# Bimodal distribution with zero inflation: CORRECT

# Difference between 0 and 4 hours, per subject----
dff <- data.table(gene = ra$gene,
                  ra[, seq(3, 13, 2), with = FALSE] - 
                    ra[, seq(2, 12, 2), with = FALSE])
colnames(dff)[-1] <- gsub(pattern = "-4h",
                          replacement = "",
                          x = colnames(dff)[-1])
dff

# Remove genes with zero change----
gene.keep <- rowSums(dff[, -1]) > 0
length(unique(dff$gene)) - sum(gene.keep)
# 14,760 genes to be removed

dff <- dff[gene.keep, ]
dff
# 11,604 genes left out of 26,364 genes left
summary(dff)

# Top differences----
rnk <- data.table(gene = dff$gene,
                  rank = rowSums(as.matrix(dff[, -1])))
rnk <- rnk[order(rnk$rank,
                 decreasing = TRUE), ]
rnk

# Separate top 20 genes with a increased expressions----
tmp <- dff[dff$gene %in% rnk$gene[1:20], ]
tmp

tmp <- melt.data.table(tmp,
                       id.vars = 1,
                       measure.vars = 2:7,
                       variable.name = "Subject_ID",
                       value.name = "Diff")
tmp <- merge(tmp,
             unique(map[, c(2, 4)]),
             by = "Subject_ID")
tmp$gene <- factor(tmp$gene,
                   levels = rev(rnk$gene[1:20]))

# Plot samples
p3 <- ggplot(tmp,
             aes(x = Diff,
                 y = gene,
                 fill = Sex,
                 group = Subject_ID)) +
  geom_point(size = 3,
             alpha = 0.8) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_x_continuous("4h - 0h TMP Difference") +
  scale_y_discrete("Gene") +
  ggtitle("Top 20 Genes with Increased Expressions")
ggplotly(p3)

# Separate top 20 genes with increased expressionsin ALL samples----
rnk2 <- data.table(gene = dff$gene,
                   pos = apply(dff[, -1],
                               1,
                               function(a) {
                                 sum(a > 0)
                               }))
table(rnk2$pos)
# 1    2    3    4    5    6 
# 989 1399 3697 3869 1437  213 

rnk <- merge(rnk, 
             rnk2,
             by = "gene")
rnk

rnk.pos <- rnk[pos == 6, ]
rnk.pos <- rnk.pos[order(rnk.pos$rank,
                     decreasing = TRUE), ]
rnk.pos
  
tmp <- dff[dff$gene %in% rnk.pos$gene[1:20], ]
tmp

tmp <- melt.data.table(tmp,
                       id.vars = 1,
                       measure.vars = 2:7,
                       variable.name = "Subject_ID",
                       value.name = "Diff")
tmp <- merge(tmp,
             unique(map[, c(2, 4)]),
             by = "Subject_ID")
tmp$gene <- factor(tmp$gene,
                   levels = rev(rnk.pos$gene[1:20]))

# Plot samples
p4 <- ggplot(tmp,
             aes(x = Diff,
                 y = gene,
                 fill = Sex,
                 group = Subject_ID)) +
  geom_point(size = 3,
             alpha = 0.8) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_x_continuous("4h - 0h TMP Difference") +
  scale_y_discrete("Gene") +
  ggtitle("Top 20 Genes with Increased Expressionsin All Samples")
ggplotly(p4)

# Separate top genes with decreased expressionsin ALL samples----
rnk3 <- data.table(gene = dff$gene,
                   neg = apply(dff[, -1],
                               1,
                               function(a) {
                                 sum(a < 0)
                               }))
table(rnk3$neg)
# 0    1    2    3    4    5 
# 1371 2225 4003 3276  712   17 
# Not a single gene had a decrease in all 6 subjects. 
rnk <- merge(rnk, 
             rnk3,
             by = "gene")
rnk

# Look at the 17 with decrease in 5 subjects----
rnk.neg <- rnk[neg == 5, ]
rnk.neg <- rnk.neg[order(rnk.neg$rank,
                         decreasing = FALSE), ]
rnk.neg

tmp <- dff[dff$gene %in% rnk.neg$gene, ]
tmp

tmp <- melt.data.table(tmp,
                       id.vars = 1,
                       measure.vars = 2:7,
                       variable.name = "Subject_ID",
                       value.name = "Diff")
tmp <- merge(tmp,
             unique(map[, c(2, 4)]),
             by = "Subject_ID")
tmp$gene <- factor(tmp$gene,
                   levels = rnk.neg$gene)

# Plot samples
p5 <- ggplot(tmp,
             aes(x = Diff,
                 y = gene,
                 fill = Sex,
                 group = Subject_ID)) +
  geom_point(size = 3,
             alpha = 0.8) +
  geom_vline(xintercept = 0, 
             linetype = "dashed") +
  scale_x_continuous("4h - 0h TMP Difference") +
  scale_y_discrete("Gene") +
  ggtitle("Top 20 Genes with Decreased Expressionsin All Samples")
ggplotly(p5)

# Look at the top genes with decrease in 4 subjects----
rnk.neg <- rnk[neg == 4, ]
rnk.neg <- rnk.neg[order(rnk.neg$rank,
                         decreasing = FALSE), ]
rnk.neg

tmp <- dff[dff$gene %in% rnk.neg$gene[1:20], ]
tmp

tmp <- melt.data.table(tmp,
                       id.vars = 1,
                       measure.vars = 2:7,
                       variable.name = "Subject_ID",
                       value.name = "Diff")
tmp <- merge(tmp,
             unique(map[, c(2, 4)]),
             by = "Subject_ID")
tmp$gene <- factor(tmp$gene,
                   levels = rnk.neg$gene[1:20])

# Plot samples
p5 <- ggplot(tmp,
             aes(x = Diff,
                 y = gene,
                 fill = Sex,
                 group = Subject_ID)) +
  geom_point(size = 3,
             alpha = 0.8) +
  geom_vline(xintercept = 0, 
             linetype = "dashed") +
  scale_x_continuous("4h - 0h TMP Difference") +
  scale_y_discrete("Gene") +
  ggtitle("Top 20 Genes with Decreased Expressionsin All Samples")
ggplotly(p5)

sessionInfo()
# sink()