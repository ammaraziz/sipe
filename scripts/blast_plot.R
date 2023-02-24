library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

files = list.files(pattern = "*all.sorted.tsv", path = "../results/blast/", full.names = T)

for (f in files) {
  tmp = read.delim(f, sep = "\t")
  counts = tmp %>% count(top_hit)
  name = gsub("_L001_R1_001.fastq.gz.blast.results", '' , f)
  plot = ggplot(counts, aes(x = fct_reorder(top_hit, n), y = n, fill = top_hit)) +
    geom_col() + 
    theme(legend.position="none", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title=name,
       x ="Serotype", y = "Reads") + 
    coord_flip()
  plot2 = ggplot(tmp, aes(y = top_hit, group = top_hit, color = top_hit, alpha = 0.1)) + 
    geom_linerange(aes(xmin = subject_start, xmax = subject_end, y = top_hit)) + 
    theme(legend.position="none")
  plot2
  
  ggsave(paste0(name, "_loc.pdf"), plot2, units = "mm", width = 210, height = 297)
  ggsave(paste0(name, "_hits.pdf"), plot, units = "mm", width = 148, height = 210)
}

# 
# tmp3 = tmp2 %>% rowwise %>% mutate(s = list(seq(V9, V10)))
# maxl = max(sapply(tmp3$s,length))
# reads_df = do.call(cbind, lapply(tmp3$s,function(x) x[1:maxl]))
# reads_df = as.data.frame(reads_df)
# names(reads_df) = paste(tmp3$V1, tmp3$V2, sep = "|")
# 
# final = reads_df %>% pivot_longer(cols = everything(),
#                           names_to = "read_name", 
#                           values_to = "positions") %>%
#   separate(col = "read_name", into = c("read_name", "blast_hit"), sep = "\\|")
# 
# read_counts = final %>% 
#   group_by(blast_hit, positions) %>% 
#   summarise(n=n()) %>%
#   filter(!is.na(positions)) %>%
#   filter(n > 10) %>%
#   separate(blast_hit, into = c("organism", "accession"), sep = "_")
# 
# ggplot(read_counts, aes(x = positions, y = n)) +
#   geom_col(aes(color = organism)) + 
#   facet_wrap(~organism)
