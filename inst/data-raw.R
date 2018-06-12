library(tximportData)
library(tximport)
library(SummarizedExperiment)
library(utils)
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples[["condition"]] <- as.factor(c(rep("c1", 3), rep("c2", 3)))
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)

# tx2gene links transcript IDs to gene IDs for summarization
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

rownames(samples) <- names(files)
metadata <- list(countsFromAbundance = txi$countsFromAbundance)

se <- SummarizedExperiment(assays = list(raw = round(txi$counts),
                                   tpm = txi$abundance,
                                   length = txi$length),
                     colData = samples,
                     metadata = metadata)
save(se, file = "data/se.rda", compress = "xz")
