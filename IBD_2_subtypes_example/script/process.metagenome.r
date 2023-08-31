library(dplyr)
library(readxl)
# read metagenome data
metagenome <- read_xlsx("IBD_2_subtypes_example/input_data/NIHMS1510763-supplement-Dataset_6.xlsx",
            sheet = 1,
            skip = 1) %>% as.data.frame()
Diagnosis <- metagenome[3, ]
sample_info <- names(Diagnosis)

#clinical info
meta <- data.frame(sample = sample_info, Diagnosis = unlist(Diagnosis))
meta <- meta[-1, ]
write.csv(meta, "IBD_2_subtypes_example/input_data/mg.meta.csv", row.names = F)


#process metagenome gene abundance data 
metagenome <- read_xlsx("IBD_2_subtypes_example/input_data/NIHMS1510763-supplement-Dataset_6.xlsx",
            sheet = 1,
            skip = 1) %>% as.data.frame()
metagenome <- metagenome[-c(1:8), ]

metagenome$`# Feature / Sample` <-
  sub(":.*", "", metagenome$`# Feature / Sample`)

#convert enzyme to keggID
ezyme <- read.table("IBD_2_subtypes_example/input_data/enzyme_ko.txt")
ezyme$ko <- sub(".*:", "", ezyme$V2)
ezyme$ec <- sub(".*:", "", ezyme$V1)
metagenome$keggid <- ezyme$ko[match(metagenome$`# Feature / Sample`, ezyme$ec)]

#
metagenome.select <- metagenome[, c(222, 2:221)]
metagenome.select[2:ncol(metagenome.select)] <-
  apply(metagenome.select[2:ncol(metagenome.select)], 2, as.numeric)
metagenome.select <-
  aggregate(metagenome.select[2:ncol(metagenome.select)],
            by = list(metagenome.select$keggid),
            FUN = sum)
row.names(metagenome.select) <- metagenome.select$Group.1
metagenome.select <- metagenome.select[, -1]

write.csv(metagenome.select, "IBD_2_subtypes_example/input_data/mg.expr.csv")

