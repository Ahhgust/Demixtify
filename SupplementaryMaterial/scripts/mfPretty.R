#!/usr/local/bin/Rsript

suppressPackageStartupMessages(require(tidyverse))

args <- commandArgs(trailingOnly=TRUE)


readFile <- function(f) {
    t <- read_tsv(f, col_names=F, col_types=cols(), progress=FALSE) %>%
        mutate(File=f)
    return(t)
}


tib <- lapply(args, readFile) %>%
    bind_rows()

colnames(tib) <- c("Label", "Fraction", "LL", "Filename")


filter(tib, Label!='nsnps') %>%
    group_by(Filename) %>%
    arrange(desc(LL), .by_group=TRUE) %>%
    summarize(MaxLL=LL[[1]],
              MF=Fraction[[1]],
              PenultLL=LL[[2]],
              PenultMF=Fraction[[2]],
              SingleSourceLL=LL[[ which(Fraction==0) ]]
              ) %>%
    ungroup() -> maxLs

filter(tib, Label=='nsnps') %>%
    select(Filename, ErrorRate=Fraction, NSnps=LL) %>%
    left_join(maxLs, by="Filename") -> maxLs

cat(format_tsv(maxLs))
