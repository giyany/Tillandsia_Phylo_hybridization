rm(list=ls())
library(plyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(data.table)
library(ape)
################################# overview #####################################

# The main data produced by Twisst is a weights file which has columns for each
# topology and their number of observations of that topology within each
# genealogy. Weights files produced by Twisst also contain initial comment lines
# speficying the topologies.

# The other data file that may be of interest is the window data. That is, the
# chromosome/scaffold and start and end positions for each of the regions or
# windows represented in the weights file.

# Both of the above files can be read into R, manipulated and plotted however
# you like, but I have written some functions to make these tasks easier.
# These functions are provided in the script plot_twisst.R

################### load helpful plotting functions #############################
setwd("/home/gil/Science/PhyloProject/Analysis/twisst/")

source("plot_twisst_mod.R")

################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.
#sorting my data file so the coordinates will make sense 

#weights file with a column for each topology
weights_file <- "output.run29.weights.csv.gz"


# It is not necessary to import window data files, but if you do there should be one for
# each weights file

#coordinates file for each window
window_data_file <- "output.run29.data.rename.tsv"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file, reorder_by_start = T)


############################## combined plots ##################################
# there are a functions available to plot both the weightings and the topologies
run="run29"
#a summary plot shows all the topologies and a bar plot of their relative weightings
plot.twisst.summary(twisst_data, lwd=4, cex=1.1, ylim=c(0,0.7))
png("run14_summary.png")

#or plot ALL the data across the chromosome(s)
#Note, this is not recommended if there are large numbers of windows.
# instead, it is recommended to first smooth the weghtings and plot the smoothed values
dev.off()


png(paste0(run,"twisst_summary_weights.png"),width = 900,height = 600,
    units = "px", pointsize = 17)
plot.twisst.summary(twisst_data, lwd=3, cex=1.2, ylim=c(0,0.9))
Sys.sleep(1)
dev.off()

######################## subset to only specific regions #########################

#regions to keep (more than one can be specified)
regions <- c("Scaffold_2199","Scaffold_658","Scaffold_2318","Scaffold_2041","Scaffold_2096","Scaffold_2315","Scaffold_2317",
             "Scaffold_2287")

regions <- c("Scaffold_326","Scaffold_1174","Scaffold_2314","Scaffold_2313","Scaffold_2225",
             "Scaffold_1383","Scaffold_819","Scaffold_15","Scaffold_2073","Scaffold_2321")

regions <- c("Scaffold_2320","Scaffold_346","Scaffold_845","Scaffold_2316",
             "Scaffold_2312","Scaffold_53","Scaffold_2319")

#subset twisst object for these
twisst_data_cont2 <- subset.twisst.by.regions(twisst_data, regions)

#or plot ALL the data across the chromosome(s)
#Note, this is not recommended if there are large numbers of windows.
# instead, it is recommended to first smooth the weghtings and plot the smoothed values
dev.off()
plot.twisst(twisst_data_cont2, mode=2, show_topos=F)

regions_small <- c("Scaffold_2199")
twisst_data_cont2 <- subset.twisst.by.regions(twisst_data, regions_small)

# make smooth weightings and plot those across chromosomes
twisst_data_smooth <- smooth.twisst(twisst_data_cont2, span_bp = 600000, spacing = 100)

plot.twisst(twisst_data_smooth, mode=2,show_topos=F,include_region_names=F) #mode 2 overlays polygons, mode 3 would stack them

prefix="run10"
# png
png(paste0(prefix,"_",regions,"_smooth_win_span_600000.png"),width = 1200,height = 2400)
plot.twisst(twisst_data_smooth, mode=3,show_topos=F,include_region_names=F) #mode 2 overlays polygons, mode 3 would stack them
Sys.sleep(1)
dev.off()

############# smoothing by windows - test Tibo


#regions to analyse
regions
listregions <- c("Chr18")
# set parameters for fitting
maximumdist=300000 # windows to consider to generate average
shiftwin=10000 # sliding approach  
minimumnumberofresults=3  # minimum number of windows to compute local average
prefix="run14"

pdf(file = paste0(prefix,"_summary_twisstplots_maxdist",maximumdist,"_shift",shiftwin,"_minwin",minimumnumberofresults,".pdf"),width = 12,height = 7)

for (region in listregions){
  print(paste0("starting to analyse ",region))

  #subset based on these regions
  twisst_data_cont2 <- subset.twisst.by.regions(twisst_data, region)
  #this can then be used in all the same plotting functions above.


  print(paste0("working on ",length(twisst_data_cont2$pos[[region]])," results from Twisst"))
  #print(head(twisst_data_cont2$pos[[region]]))#Scaffold_1174)
  #print(nrow(twisst_data_cont2$weights[[region]]))
  weightstransf=NULL
  for (posindex in 1:length(twisst_data_cont2$pos[[region]])){
    weightstransf=rbind(weightstransf,cbind(twisst_data_cont2$pos[[region]][posindex],twisst_data_cont2$weights[[region]][posindex,]))
  }
  weightstransf=as.data.frame.matrix(weightstransf)
  colnames(weightstransf)<- c("pos",colnames(twisst_data_cont2$weights[[region]]))
  
  # loop
  dataweightstransflocal <- data.frame(matrix(ncol = ncol(weightstransf)+2, nrow = 0))
  #dataweightstransflocal=NULL
  colnames(dataweightstransflocal)<-c("startwin","endwin",colnames(weightstransf))
  #print(colnames(dataweightstransflocal))
  for (windows in floor(maximumdist/shiftwin):ceiling(max(weightstransf$pos)/shiftwin)){ # for all windows from the first positive startwin values
    #print(windows)
    startwin=shiftwin*windows-maximumdist+1
    endwin=shiftwin*windows
    weightstransflocal<-subset(weightstransf, (weightstransf$pos >= startwin) & (weightstransf$pos <= endwin ))
    if (nrow(weightstransflocal) >= minimumnumberofresults ){
      averageweights=as.data.frame.matrix(cbind(startwin,endwin,t(apply(weightstransflocal,MARGIN=2,FUN=mean))))
      colnames(averageweights)<-c("startwin","endwin",colnames(weightstransf))
      #print(averageweights)
      dataweightstransflocal<-rbind.fill(dataweightstransflocal,averageweights)
      #print(cbind(startwin,endwin,averageweights))
    }
    else {
      averageweights=as.data.frame.matrix(cbind(startwin,endwin, t(rep("NA",ncol(weightstransf)))))
      colnames(averageweights)<-c("startwin","endwin",colnames(weightstransf))
      #print(averageweights)
    dataweightstransflocal<-rbind(dataweightstransflocal,averageweights)
    }
  }
  # output the matrix
  write.table(x=dataweightstransflocal,file = paste0(prefix,"_",region,"_avgweightings_win_",maximumdist,"_slidwin_",shiftwin,".txt"),row.names = FALSE, quote = FALSE,sep = "\t")

  # transform the matrix to produce a plot
  newdataweightstransflocal=NULL
  for (column in 4:ncol(dataweightstransflocal)){
    newdataweightstransflocal=as.data.frame.matrix(rbind(newdataweightstransflocal,cbind(dataweightstransflocal[,1],dataweightstransflocal[,2],colnames(dataweightstransflocal[column]),noquote(dataweightstransflocal[,column]))))
  }
  colnames(newdataweightstransflocal) <- c("startwin","endwin","topology","weights")
  write.table(x=newdataweightstransflocal,file = paste0(prefix,"_",region,"_avgweightings_win_",maximumdist,"_slidwin_",shiftwin,"prggplot2.txt"),row.names = FALSE, quote = FALSE,sep = "\t")

  data2plot <- read.table(file = paste0(prefix,"_",region,"_avgweightings_win_",maximumdist,"_slidwin_",shiftwin,"prggplot2.txt"),header = TRUE, sep = "\t")

  #print(data2plot)
  missingdata2plot <- subset(data2plot,data2plot$topology=="topo1" & is.na(data2plot$weights==TRUE))
  toplot <- ggplot(data2plot, aes(startwin, weights, colour=topology)) +  geom_bar(width=0.8, stat="identity") + 
    scale_colour_manual(na.value = "black", values = c(
      "#0075DC", #Blue
      "#2BCE48", #Green
      "#FFA405", #Orpiment
      "#5EF1F2", #Sky
      "#FF5005", #Zinnia
      "#005C31", #Forest
      "#00998F", #Turquoise
      "#FF0010", #Red
      "#9DCC00", #Lime
      "#003380", #Navy
      "#F0A3FF", #Amethyst
      "#740AFF", #Violet
      "#426600", #Quagmire
      "#C20088", #Mallow
      "#94FFB5")) +#Jade
    xlab("position")+ ylab("Averaged Weights")+theme(legend.position="none")+ggtitle(region)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
          axis.text.x = element_text(colour="black",size=25,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(colour="black",size=25,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.2,face="italic"),
          axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="italic"),
          plot.title = element_text(size=30,hjust = 0.5)) +
  #increase number of ticks on x axis
  scale_x_continuous(breaks = round(seq(min(0), max(data2plot$startwin), by = 5000000),0)) 
    # if the "scale_colour_manual(na.value = "black"" does not work well the line below is needed 
#    toplot <- toplot + geom_segment(data=missingdata2plot,aes(x=startwin,y = 0, xend=startwin,yend = 1), colour = "black") 
    # print the plot
  # file for plots
    print(toplot)
    # png
    png(paste0(prefix,"_",region,"_avgweightings_win_",maximumdist,"_slidwin_",shiftwin,"_minwin",minimumnumberofresults,".png"),width = 1200,height = 600,
        units = "px", pointsize = 17)
    print(toplot)
    Sys.sleep(1)
    dev.off()
}
dev.off()


##################### individual plots: raw weights ############################

#plot raw data in "stepped" style, with polygons stacked.
#specify stepped style by providing a matrix of starts and ends for positions
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot.weights(weights_dataframe=twisst_data$weights[[1]], positions=twisst_data$window_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)

#plot raw data in stepped style, with polygons unstacked (stacked =FLASE)
#use semi-transparent colours for fill
plot.weights(weights_dataframe=twisst_data_contig0$weights[[1]], positions=twisst_data_contig0$window_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=paste0(topo_cols,80), stacked=FALSE)


#################### individual plots: smoothed weights ########################

#plot smoothed data with polygons stacked
plot.weights(weights_dataframe=twisst_data_smooth$weights[[1]], positions=twisst_data_smooth$pos[[1]],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)

#plot smoothed data with polygons unstacked
plot.weights(weights_dataframe=twisst_data_smooth$weights[[1]], positions=twisst_data_smooth$pos[[1]],
             line_cols=topo_cols, fill_cols=paste0(topo_cols,80), stacked=FALSE)


########################### plot topologies using Ape ##########################
#unrooted trees
for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- ladderize(unroot(twisst_data$topos[[i]]))

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "unrooted", edge.color=topo_cols[n], edge.width=5, rotate.tree = 90, cex = 1, adj = .5, label.offset=.2)
  mtext(side=3,text=paste0("topo",n))
  }


#rooted topologies
for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- root(twisst_data$topos[[i]], "D", resolve.root = T)

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "clad", edge.color=topo_cols[n], edge.width=5, label.offset=.1, cex = 1)
  mtext(side=3,text=paste0("topo",n))
  }

#################### subset to only the most abundant topologies #################

#get list of the most abundant topologies (top 5 in this case)
top8_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[1:8]

#subset twisst object for these
twisst_data_top8topos <- subset.twisst.by.topos(twisst_data, top8_topos)
#this can then be used in all the same plotting functions above.

contig2_top5<-subset.twisst.by.regions(twisst_data_top5topos,regions)

######################## subset to only specific regions #########################

#regions to keep (more than one can be specified)
regions <- c("Scaffold_1174")

#subset twisst object for these
twisst_data_cont2 <- subset.twisst.by.regions(twisst_data, regions)
#this can then be used in all the same plotting functions above.
#regions to keep (more than one can be specified)
regions <- c("Scaffold_2199","Scaffold_658","Scaffold_2318","Scaffold_2041")

#subset twisst object for these
twisst_data_contig0 <- subset.twisst.by.regions(twisst_data, regions)

######################### my way ###################################

#my way more simple twisst plots
library(dplyr)
library(ggplot2)
library(data.table)

weights_data<-read.table("output.run7.weights.csv.gz",head=T)
window_data <- read.table("output.run7.data.tsv",head=T)

#calc weight for each topo
weights_data$top1_weight<-weights_data$topo1/rowSums(weights_data[1:15])
weights_data$top2_weight<-weights_data$topo2/rowSums(weights_data[1:15])
weights_data$top3_weight<-weights_data$topo3/rowSums(weights_data[1:15])
weights_data$top4_weight<-weights_data$topo4/rowSums(weights_data[1:15])
weights_data$top5_weight<-weights_data$topo5/rowSums(weights_data[1:15])
weights_data$top6_weight<-weights_data$topo6/rowSums(weights_data[1:15])
weights_data$top7_weight<-weights_data$topo7/rowSums(weights_data[1:15])
weights_data$top8_weight<-weights_data$topo8/rowSums(weights_data[1:15])
weights_data$top9_weight<-weights_data$topo9/rowSums(weights_data[1:15])
weights_data$top10_weight<-weights_data$topo10/rowSums(weights_data[1:15])
weights_data$top11_weight<-weights_data$topo11/rowSums(weights_data[1:15])
weights_data$top12_weight<-weights_data$topo12/rowSums(weights_data[1:15])
weights_data$top13_weight<-weights_data$topo13/rowSums(weights_data[1:15])
weights_data$top14_weight<-weights_data$topo14/rowSums(weights_data[1:15])
weights_data$top15_weight<-weights_data$topo15/rowSums(weights_data[1:15])

my_twisst_data=cbind(window_data,weights_data)

my_twisst_data_melted<-melt(as.data.table(my_twisst_data),id.vars=c("scaffold","start","end"), measure.var=patterns("_weight"),value.name = "weight", variable.name="topo")

my_color_pal<-c("#0075DC", "#2BCE48", "#FFA405", "#5EF1F2", "#FF5005", "#005C31",
"#00998F", #Turquoise
"#FF0010", #Red
"#9DCC00", #Lime
"#003380", #Navy
"#F0A3FF", #Amethyst
"#740AFF", #Violet
"#426600", #Quagmire
"#C20088", #Mallow
"#94FFB5") #Jade

#just scaffold 2199
twisst_data_melted_scaffold_2199<-subset(my_twisst_data_melted, scaffold=="Scaffold_2199")
twisst_data_melted_scaffold_658<-subset(my_twisst_data_melted, scaffold=="Scaffold_658")
twisst_data_melted_scaffold_2318<-subset(my_twisst_data_melted, scaffold=="Scaffold_2318")
twisst_data_melted_scaffold_2041<-subset(my_twisst_data_melted, scaffold=="Scaffold_2041")

ggplot(twisst_data_melted_scaffold_2199, aes(fill=topo, y=weight, x=start)) + 
  geom_bar(position="stack", stat="identity",width=30000) + scale_fill_manual(values=my_color_pal) + theme_minimal()

p2<- ggplot(twisst_data_melted_scaffold_658, aes(fill=topo, y=weight, x=start)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=my_color_pal) + theme_minimal()

p3<- ggplot(twisst_data_melted_scaffold_2318, aes(fill=topo, y=weight, x=start)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=my_color_pal) + theme_minimal()
p4<- ggplot(twisst_data_melted_scaffold_2041, aes(fill=topo, y=weight, x=start)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=my_color_pal) + theme_minimal()
plot_grid(p1,p2, nrow = 2)


p1
