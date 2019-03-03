#!/usr/local/bin/Rscript

library(plyr)
library(ggplot2)
library(reshape2)
library(grid)	
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
output_directory = args[1]
phase_cmdline = args[2]
cycle = args[3]
phasing_score_files = list.files(path = output_directory, pattern = "*.phasing_score")
abundance_files = list.files(path = output_directory, pattern = "*.abundance")


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) 
{
	plots <- c(list(...), plotlist)
	numPlots = length(plots)
	if (is.null(layout)) 
	{
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	if (numPlots==1) 
	{
		print(plots[[1]])
	} 
	else 
	{
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		for (i in 1:numPlots) 
		{
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

graph_list<-list()
i<-1


print("This are the phasing score files")
print(phasing_score_files)
phasing_score_files<-phasing_score_files[1:length(phasing_score_files)]

for(file_num in seq(1:length(phasing_score_files)))
{
	phasing_score_file = phasing_score_files[file_num]
	abundance_file = abundance_files[file_num]
	
	splitted_string <- strsplit(strsplit(abundance_file,split = ".abundance")[[1]],split="_")[[1]]
	phase=as.numeric(splitted_string[1])
	chromosome=splitted_string[2]
	start=as.numeric(splitted_string[3])
	end=as.numeric(splitted_string[4])
	phasing_score_file = paste0(output_directory,"/",phasing_score_files[file_num])
	abundance_file = paste0(output_directory,"/",abundance_files[file_num])
	colour_string<-paste(chromosome,"->",toString(start),"-",toString(end))
	
	if(is.na(start) && is.na(end))
		next
	if(start==0 && end==0)
		next
	print(phasing_score_file)
	print(abundance_file)
	print(colour_string)
	
	
	score<-read.table(phasing_score_file,head=F,sep='\t')
	names(score)<-c('coordinate','score');
	abundance<-read.table(abundance_file,head=F,sep='\t')
	names(abundance)<-c('coordinate','abundance');
	
	p1<-ggplot(score,aes(x=coordinate,y=score))+
			geom_line()+
			scale_x_continuous(limits=c(start,end),breaks=seq(start,end,phase))+
			theme(axis.text.x=element_text(angle=90,size=12,vjust=0.5),axis.text.y=element_text(size=15))+
			ggtitle("Phasing score",subtitle=colour_string)+
			theme(plot.title = element_text(hjust=0.5,face="bold"),
					plot.subtitle = element_text(hjust=0.5));
	
	p2<-ggplot(abundance,aes(x=coordinate,y=abundance))+
			geom_point(size=1)+scale_x_continuous(limits=c(start,end),breaks=seq(start,end,phase))+
			theme(axis.text.x=element_text(angle=90,size=12,vjust=0.5),axis.text.y=element_text(size=15))+
			ggtitle("Abundance",subtitle=colour_string)+
			theme(plot.title = element_text(hjust=0.5,face="bold"),
					plot.subtitle = element_text(hjust=0.5));
	
	graph_list[[i]]<-p1
	graph_list[[i+1]]<-p2
	#multiplot(p1,p2)
	#dev.off()
	i<-i+2
}

multi_plot <- marrangeGrob(grobs=graph_list, nrow=1, ncol=1,top=NULL)
ggsave(paste0(output_directory,"/","plots","_phase_",phase_cmdline,"_cycle_",cycle,".pdf"), multi_plot, width=8.5)


