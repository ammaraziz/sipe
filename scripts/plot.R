library(ggplot2)
files = list.files(path = "../analysis/depth/raw/", pattern = ".txt",full.names = T)

for (f in files) {
  output_name = gsub(".txt", ".pdf", f)
  
  dat = read.delim(f, header = F, sep = "\t")
  names(dat) = c("genome", "pos", "depth")
  p = ggplot(dat, aes(x = pos, y = depth, group = genome, color = genome)) + 
    geom_col() + facet_wrap(~genome) +
    theme(legend.position="bottom", legend.box = "horizontal")
  ggsave(output_name, p, device = pdf, width = 210, height = 297, units = 'mm') 
}
