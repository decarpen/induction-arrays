ggheat = function(dat, colors=c("blue","white","red"), dist.method='spearman',
				  clust.method='average', file.name=NULL, width.adjust=2.1) {
  #Clusters and makes a heatmap of a matrix and saves it to a PDF.
  #
  #Args:
  #  dat: a numeric matrix of your data. Row and column names will be used as labels.
  #  colors: a vector of colors to use, the first value is low, second is the midpoint, and third is high.
  #  dist.method: distance metric calculation. Can be 'euclidean' or 'spearman'
  #  clust.method: 'complete' or 'average'
  #  file.name: The name for the saved plot. (E.g., "heatmap.pdf")
  #  width.adjust: 
  
  #Returns:
  #  g: a ggplot instance of the heatmap
  
  require(ggplot2)
  require(cba)
  
  	if (dist.method=='euclidean') {
  		d = dist(dat)
  	}
  	else {
  		d = 1-cor(dat,method='spearman')
  		d = as.dist(d)
  	}
  	
  	if (clust.method == 'complete') {
  		dat = dat[hclust(d,method='complete')$order,]
  	}
  	else {
  		dat = dat[hclust(d,method='average')$order,]
  	}

  	#rownames(dat) = hclust(dist(dat))$order
	rows = dim(dat)[1]
	cols = dim(dat)[2]
	dat.m = cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols,each=rows), melt(dat))
  

	win.width=(ncol(dat)/8)+width.adjust
	win.height=(nrow(dat)/12)

	base_size=8

	g = ggplot(data=dat.m) + 
		geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd,fill=value,colour=value))
	
	g = g+scale_x_continuous(breaks=(1:cols)-0.5,labels=colnames(dat)) +
		scale_y_continuous(breaks=(1:rows)-0.5,labels=rownames(dat)) +
		opts(panel.grid.minor=theme_line(colour=NA), panel.grid.major=theme_line(colour=NA),
			panel.background=theme_rect(fill="light grey", colour=NA), axis.text.y=theme_text(size=6, hjust=1),
			axis.ticks=theme_blank(),axis.text.x=theme_text(angle=90, hjust=1,size=6), axis.title.x=theme_blank()) +
        scale_fill_gradient2(low=colors[1],mid=colors[2],high=colors[3],midpoint=0, name="Fold Change") +
  		scale_colour_gradient2(low=colors[1],mid=colors[2],high=colors[3],midpoint=0, name="Fold Change")

  if (!is.null(file.name) {
  	ggsave(filename=file.name, width=win.width, height=win.height)
  }
  
  return(g)
}