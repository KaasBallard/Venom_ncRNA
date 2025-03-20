#set environment
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Scripts/R/x_misc/HierachicalEdgeModeling')

#open files and create the df objects to generate the network
edges1 <- read.csv("./aca_E_ggraph.txt", header=T, sep="\t", as.is=T)
connect1 <- read.csv("./aca_CONNECT_ggraph.txt", header=T, sep="\t", as.is=T) %>% 
  filter(str_detect(from, 'mir|let')) %>% 
  filter(!str_detect(to, 'mir|let'))

#create the vertices1 object
vertices1 <- data.frame(name = unique(c(as.character(edges1$from), as.character(edges1$to))), value = rep(1, length(unique(c(as.character(edges1$from), as.character(edges1$to))))))
vertices1$group  <-  edges1$from[ match( vertices1$name, edges1$to ) ]

#Create the labels and set the angles
vertices1 <- vertices1[order(vertices1$group),]
vertices1$id <- NA
myleaves <- which(is.na( match(vertices1$name, edges1$from) ))
nleaves <- length(myleaves)
vertices1$id[ myleaves ] <- seq(1:nleaves)
vertices1$angle <- 90 - 360 * vertices1$id / nleaves
#calculate the alignment of labels: right or left; If I am on the left part of the plot, my labels have currently an angle < -90
vertices1$hjust <- ifelse( vertices1$angle < -90, 1, 0)
#flip angle BY to make them readable
vertices1$angle <- ifelse(vertices1$angle < -90, vertices1$angle+180, vertices1$angle)

#create the igraph object and the to and from objects
mygraph <- igraph::graph_from_data_frame( edges1, vertices=vertices1 )
from  <-  match( connect1$from, vertices1$name)
to  <-  match( connect1$to, vertices1$name)


#plotting network
plot <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.25, width=0.5, aes(colour=after_stat(index))) +
  scale_edge_colour_distiller(palette = "RdPu") +
  
  geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust, colour=group), size=1.5, alpha=1) +
  
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2)) +
  scale_colour_manual(values= c("darkgoldenrod4","grey14","orange","aquamarine4")) +
  scale_size_continuous( range = c(0.1,2) ) +
  
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
plot
#END






