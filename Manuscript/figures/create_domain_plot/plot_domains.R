### SJS 10/22/14 ###
### R script to plot figure domains_naive_struc.pdf. Requires ggplot2 and cowplot (see here - https://github.com/wilkelab/cowplot )
### Uses data files created in python script prepare_domain_plot.py


library(cowplot) # requires ggplot so automatically loaded here.

# Set up domain colors.
palette <- c("-" = "gray85", "A" = "white", "O" = "navy", "I" = "dodgerblue1", "M" = "firebrick2") #I:"#A3D5FB", #O:"#0095FF"

# Read in naive and structural aln, consensus data frames.
aln_n <- read.table('naive_plotaln.txt', header=T)
cons_n  <- read.table('naive_plotcons.txt', header=T)
aln_s <- read.table('struc_plotaln.txt', header=T)
cons_s  <- read.table('struc_plotcons.txt', header=T)

# Create new data frames for use with geom_rect.
df_aln_n=data.frame(x1=aln_n$x, x2=aln_n$x + 20, y1=aln_n$y, y2=aln_n$y + 20, t=aln_n$domain)
df_cons_n=data.frame(x1=cons_n$x, x2=cons_n$x + 20, y1=cons_n$y - 12, y2=cons_n$y - 2, t=cons_n$domain)
df_aln_s=data.frame(x1=aln_s$x + 1970, x2=aln_s$x + 1990, y1=aln_s$y, y2=aln_s$y + 20, t=aln_s$domain)
df_cons_s=data.frame(x1=cons_s$x + 1970, x2=cons_s$x + 1990, y1=cons_s$y - 12, y2=cons_s$y - 2, t=cons_s$domain)

# Massive ggplot object
p <- ggplot() + geom_rect(df_aln_n, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = t), linetype=0) + 
geom_rect(df_cons_n, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = t), linetype=0) + 
geom_rect(df_aln_s, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = t), linetype=0) + 
geom_rect(df_cons_s, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = t), linetype=0) + 
scale_fill_manual(values=palette) + theme_nothing() + 
draw_text("Naive MSA", x = 880,  y = 153, size = 25, fontface = "bold") + 
draw_text("Curated MSA", x = 2420, y = 153, size = 25, fontface = "bold") +
geom_segment(aes(x = -30, y = 0, xend = -30, yend = 148), size = 1.5) +
draw_text("1", x  = -70, y = 5, size = 22, fontface = "bold") + 
draw_text("50", x  = -70, y = 55, size = 22, fontface = "bold") + 
draw_text("100", x  = -80, y = 105, size = 22, fontface = "bold") + 
geom_segment(aes(x = 1940, y = 0, xend = 1940, yend = 148), size = 1.5) +
draw_text("1", x  = 1900, y = 5, size = 22, fontface = "bold") + 
draw_text("50", x  = 1900, y = 55, size = 22, fontface = "bold") + 
draw_text("100", x  = 1890, y = 105, size = 22, fontface = "bold") + 

draw_text("Sequence", angle = 90, x = -160, y = 75, size = 25, fontface = "bold") +
draw_text("Site", x = 900,  y = -20, size = 25, fontface = "bold") + 
draw_text("Site", x = 2440, y = -20, size = 25, fontface = "bold") +
draw_text("Consensus", x = 2980, y = -7, size = 25, fontface = "bold") + 
draw_plot_label(c("A", "B"), c(-130, 1850), c(160, 160), size = 30)

# Save!
ggsave('../domains_naive_struc.png', p , width=25, height=10)

