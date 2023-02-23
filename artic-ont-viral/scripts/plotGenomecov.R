#!/usr/bin/env Rscript

# install and/or load packages
repos='http://cran.us.r-project.org'

if(!require("dplyr")){
  install.packages('dplyr', repos = repos)
  library(dplyr)
};

if(!require("ggplot2")){
  install.packages('ggplot2', repos = repos)
  library(ggplot2)
};

if(!require("scales")){
  install.packages('scales', repos = repos)
  library(scales)
};


# define command line arguments
args = commandArgs(trailingOnly=TRUE)

# check for arguments
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  input = args[1]
}

# read coverage file
genome.cov <- read.table(input, header=FALSE, 
                         sep="\t", na.strings="NA", 
                         dec=".", strip.white=TRUE)

# rename the column names (header)
genome.cov <- genome.cov %>% dplyr::rename(Chr="V1", Pos="V2", depth="V3")


# define output name
sample <- strsplit(basename(input), ".", fixed=T)[[1]][1]

# plot
pdf(paste(sample, '.genome.coverage.pdf', sep = ''), width=15.25, height=9.5)
cov.plot <- ggplot(data = genome.cov, aes(x=Pos, y=depth)) +
  geom_ribbon(aes(ymin=0, ymax=depth), fill="#9ACD32", data=) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=10^c(0:10),
                     labels=trans_format('log10', math_format(10^.x)),
                     expand=c(0, 0)) +
  
  expand_limits(y=1) +
  ylab(bquote('log'[10]~'(Coverage+1)')) +
  xlab('Position (bp)') +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    panel.background = element_rect(colour = "black", fill=NA, size=1.0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.0))
print(cov.plot)
dev.off()


# plot(genome.cov$Pos, log10(genome.cov$depth), type = "l", 
#      xlab = "Position (bp)", 
#      ylab = bquote('log'[10]~'(Coverage)'))
# polygon(log10(genome.cov$depth), col = "steelblue")
# 
# 
# library(lattice)
# 
# xyplot(depth ~ Pos, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=genome.cov, main="depth by locus - Chr2 (SampleMA605)")
# 
# 
# d <- density(genome.cov$depth)
# plot(d, main="Kernel Density of Miles Per Gallon")
# polygon(d, col="red", border="blue")