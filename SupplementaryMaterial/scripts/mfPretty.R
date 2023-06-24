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


# mixture hypothesis
filter(tib, Label!='nsnps') %>%
    filter(Fraction>0) %>%
    group_by(Filename) %>%
    arrange(desc(LL), .by_group=TRUE) %>%
    summarize(Mixture_MaxLL=LL[[1]],
              MF=Fraction[[1]]) %>%
    ungroup() -> maxLs

# single source hypothesis
filter(tib, Label!='nsnps') %>%
    filter(near(Fraction,0)) %>%  # this should be one row/file.
    group_by(Filename) %>%
    arrange(desc(LL), .by_group=TRUE) %>%
    summarize(SingleSource_LL=LL[[1]]) %>%
    ungroup() -> maxLsNull

# and grab the error rate/number of SNPs
filter(tib, Label=='nsnps') %>%
    select(Filename, ErrorRate=Fraction, NSnps=LL) %>% # and join it all back together!
    left_join(maxLs, by="Filename") %>%
    left_join(maxLsNull, by='Filename') %>%
    mutate(LLR=Mixture_MaxLL-SingleSource_LL)  -> maxLs
#print to stdout
cat(format_tsv(maxLs))

# and make a (not so pretty) plot
filter(tib, Label!='nsnps') %>%
    ggplot(aes(x=Fraction, y=LL)) +
    geom_point() +
    facet_wrap(~Filename, ncol=1) +
    theme_bw(base_size=20) -> pl

ggsave("plotty.png", pl, height=3 + 10*length(args), width=12)

