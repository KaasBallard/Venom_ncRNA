library(tidyverse)
library(circlize)
library(rlist)
library(officer)
library(rvg)

setwd("/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Scripts/R/circos_simple_example/")

## load scaffold sizes
# will be used to initiate the sizes of each sector of the
# circos plot
scaffsize = read.delim('scaffold_sizes.txt', header = F)
# rename columns in scaffsize
names(scaffsize) = c('Chrom','size','genome_position')
# create a column of zeros, to indicate the starting
# value of each sector
scaffsize = scaffsize %>% mutate(start = 0)

## load recombination map
recomb = read.delim('viridis.recomb.bpen10.windowed.100kb.centromereMask.txt') %>%
    rename(Chrom = chrom)

## load LD data
ld = read.delim('LD_example.txt', sep = ' ')

circos.clear()
## set basic parameters for the circos plot
circos.par('track.height' = 0.2,
	cell.padding = c(0.01, 0, 0.01, 0),
	start.degree = 90,
	gap.degree = c(rep(1, 17), 20))
## initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffsize$Chrom, xlim = cbind(scaffsize$start, scaffsize$size))

## plot recombination
# here we will plot a track, each call of "circos.track" creates a new track
# "track.height" specify the size of the track
# ylim specify the vertical limits of the track, here I'm using "range" on the data I will plot to get the minimum and
# maximum value directly

# circos.track works like a loop and calls the function specified in panel.fun for each sector
# to access the current sector parameters, use CELL_META, for example to get the name of the current sector
# (i.e. the name of the chromosome) use CELL_META$sector.index, to get the y axis limit: CELL_META$ylim
# a lot of sector parameters can be accessed that way.

# inside the panel.fun function, you can call circos functions to plot stuff, but also filter your data in
# different ways depending on the sector

circos.track(ylim = range(recomb$mean, na.rm = T), track.height = 0.20,
	panel.fun = function(x, y) {
		# filter data for current chromosome
		r = recomb %>% filter(Chrom == CELL_META$sector.index)
		
		# plot genome wide recombination rate
		circos.lines(r$end, r$mean,
			col = 'black', pch = 16, cex = 0.3, area = F, baseline = 0, border = NA)
		# plot chromosome names
		# mm_y() can be used to nudge something on the y axis
		circos.text(CELL_META$xcenter, 
			CELL_META$cell.ylim[2] + mm_y(15), 
			gsub('scaffold-','',CELL_META$sector.index),
			facing = 'downward')
})


# Plot each LD one by one with a for loop

for(i in 1:nrow(ld))
{
	# need to specify sector1, position1 and sector2, position2 to plot the LD line
	circos.link(ld$CHR1[i], ld$POS1[i], ld$CHR2[i], ld$POS2[i], col = add_transparency(ld$col[i], 0.95))
}
