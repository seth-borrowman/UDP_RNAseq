library(tidyverse)

gtex <- read_delim("gtex_reads.gct",
                   delim = '\t', col_names = T, skip = 2) %>%
    rename(gene_name = Name)


udp714051 <- read_delim("quants/714051_quant/quant.sf",
                        delim = '\t', col_names = T) %>%
    select(Name, TPM) %>%
    rename(TPM = udp714051) %>%
    rename(Name = gene_name)

gtex2 <- gtex %>% 
    rowwise() %>%
    mutate(mean_tpm = mean(c_across(3:757))) %>%
    mutate(sd_tpm = sd(c_across(3:757)))

gtex_means <- gtex2 %>% select(c(Name, gene_name, mean_tpm, sd_tpm)) %>%
    left_join(., udp714051, by="gene_name") %>%
    na.omit() %>%
    mutate(udp714051_delta = (udp714051 - mean_tpm)) %>% #calculate deltas
    ungroup() %>%
    filter(mean_tpm > 0.1) %>%
    filter(udp714051 > 0.1)

sig_714051 <- gtex_means %>% filter(abs(udp714051_delta) > 1500)

scatter_714051 <- ggplot2::ggplot(data = gtex_means, aes(x=mean_tpm, y=udp714051, label = gene_name)) + 
    geom_point(size=0.6, color="darkgray", alpha=0.8) + 
    geom_point(data = sig_714051, aes(x=mean_tpm, y=udp714051), size=0.8, color="red", alpha=.9 ) +
    geom_smooth(method = "lm") + 
    scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') +
    geom_text_repel(data = sig_714051, 
                    max.overlaps = 30, 
                    size = 2.5, 
                    segment.size=0.1, 
                    nudge_x=0.06, 
                    nudge_y=0.06,
                    max.iter= 200,
                    point.padding = 0.15, 
                    segment.alpha = 1, 
                    box.padding=.15,
                    min.segment.length = unit(0.15, 'lines'),size=2.5)

png(file = "scatterplot.png", width = 1080, height = 1080)
scatter_714051
dev.off()

barplot_714051 <- ggbarplot(sig_714051, x='gene_name', y='pb714051_delta', 
                            sort.val='asc', 
                            sort.by.groups = FALSE, 
                            x.text.angle=90,
                            ylab = "Patient 714051 Delta From GTEX Mean (TPM)",
                            xlab = "Gene ID", 
                            lab.size = 1,
                            rotate = TRUE,
                            ggtheme = theme_minimal())

png(file = "barplot", width = 1080, height = 1080)
barplot_714051
dev.off()