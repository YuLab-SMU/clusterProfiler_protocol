

library(readxl)
library(dplyr)
#clinical info
metabolism <-
  read_xlsx(
    "../input_data/NIHMS1510763-supplement-Dataset_2.xlsx",
    sheet = 1,
    skip = 1
  ) %>% as.data.frame()
#clinical info
meta <- data.frame(sample = c(names(metabolism)),
             Diagnosis = c(metabolism[2, ]) %>% unlist())
meta <- meta[-1, ]
write.csv(meta, "../input_data/metabolism_meta.csv", row.names = F)

#metabolism abundance data

meta.b <- read_excel("../input_data/NIHMS1510763-supplement-Dataset_1.xlsx",
                    sheet = 1,
                    skip = 1) %>% as.data.frame()
meta.b <- meta.b[, c(1, 6)]
metabolism$`# Feature / Sample` <-
  meta.b$`Exact Match to Standard (* = isomer family)`[match(metabolism$`# Feature / Sample`, meta.b$`Metabolomic Feature`)]

#conver metabolites to KEGID
convert.name <- read.csv("../input_data/name_map.csv")
metabolism$keggid <- convert.name$KEGG[match(metabolism$`# Feature / Sample`, convert.name$Query)]

metabolism <- metabolism[!(metabolism$`# Feature / Sample` %>% is.na()), ]
metabolism <- metabolism[-c(1:7), ]

metabolism.select <- metabolism[, c(222, 2:221)]
metabolism.select[2:ncol(metabolism.select)] <-
  apply(metabolism.select[2:ncol(metabolism.select)], 2, as.numeric)
metabolism.select <-
  aggregate(metabolism.select[2:ncol(metabolism.select)],
            by = list(metabolism.select$keggid),
            FUN = sum) %>% na.omit()
row.names(metabolism.select) <- metabolism.select$Group.1
metabolism.select <- metabolism.select[, -1]

write.csv(metabolism.select, "../input_data/metabolism_expr.csv")
