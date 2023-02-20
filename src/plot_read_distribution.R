#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(extrafont)
library(ggplot2)
library(scales)

summary_file <- snakemake@input[["seqsum"]]
plot_file <- snakemake@output[["plot"]]

# read the data
names(summary_file) <- basename(dirname(summary_file))
bc_results_list <- lapply(summary_file, fread)
bc_results <- rbindlist(bc_results_list, idcol = "fc", fill = TRUE)

# release memory
rm(bc_results_list)
gc()

# draw plot
gp <- ggplot(bc_results, 
             aes(x = sequence_length_template,
                 y = mean_qscore_template)) +
  facet_wrap(~ fc) +
  theme_minimal(base_family = "Lato") +
  ylab("Mean Q score") + xlab("Read length") +
  scale_fill_viridis_c(
    # trans = log_trans(),
    # breaks = log_breaks(),
    guide = guide_colourbar(
      title = "Number of reads"
    )) +
  scale_x_continuous(trans = log_trans(base = 4),
                     breaks = trans_breaks(function(x) log(x, 4),
                                           function(x) 4^x)) +
  geom_hex() +
  geom_hline(yintercept = 9,  linetype = 2)

ggsave(plot_file,
       gp,
       width = 10,
       height = 7.5,
       units = "in",
       device=cairo_pdf)

# log
sessionInfo()