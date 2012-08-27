zeroNorm <- function(x) {
  #Zero-normalizes a matrix of array data.
  #
  #Args:
  #  x: a matrix of array data arranged with columns in increasing time order.
  #     The first column should be the time zero column.
  #
  #Returns:
  #  y: a matrix of zero-normalized data.
  
  x.zero = x[,1]
  
  y = apply(x, 2, function(a,b) {return(a-b)}, x.zero)
  
  return(y)
  
}

#--------------------------------------------------------------------------------
dropLow <- function(x, min) {
  #Removes all data that changes less than the min value.
  #
  #Args:
  #  x: the data frame to be filtered
  #  min: the minimum change threshold
  #
  #Returns:
  #  y: filtered matrix
  
  #find columns which are text
  sel.text = sapply(x[1,], is.numeric)
  
  text = x[,!sel.text]
  nums = x[,sel.text]
  
  max.list = apply(abs(nums),1,max,na.rm=T)
  sel = sapply(max.list, function(x) {x > min})
  return(x[sel,])
}

#--------------------------------------------------------------------------------
fixNames <- function(x) {
	#Cleans up annotations downloaded from PUMAdb. (Turns "YDR477W || SNF1" into "YDR477W")
	#
	#Args:
	#  x: a vector of character strings to be cleaned.
	#
	#Returns:
	#  orf: the vector of cleaned strings.
	
  orf=NULL
  x.split = strsplit(as.character(x)," ||",fixed=TRUE)
  for (i in 1:length(x.split)) {
    orf[i]=x.split[[i]][1]
  }
  return(orf)
}

#--------------------------------------------------------------------------------
orf2name <- function(x) {
	#Uses Bioconductor database "org.Sc.sgd.db" to convert ORF names to standard gene names.
	#
	#Args:
	#  x: a vector of ORF names
	#
	#Returns:
	# a list of standard gene names and ORF names (for those without standard names)
  
  require(org.Sc.sgd.db)
  
  d = org.Sc.sgdGENENAME
  mapped_genes <- mappedkeys(d)
  dd = as.list(d[mapped_genes])
  return(sapply(x, function(x) {
    if(is.null(dd[as.character(x)][[1]])) {
      return(x)
    }else
    return(dd[as.character(x)][[1]])
    }))
}

#--------------------------------------------------------------------------------
orf2desc <- function(x) {
  	#Uses Bioconductor database "org.Sc.sgd.db" to convert ORF names to gene descriptions.
	#
	#Args:
	#  x: a vector of ORF names
	#
	#Returns:
	# a list of gene descriptions for those names
  
  require(org.Sc.sgd.db)
  
  d = org.Sc.sgdDESCRIPTION
  mapped_probes <- mappedkeys(d)
  dd = as.list(d[mapped_probes])
  return(sapply(x, function(x) {dd[as.character(x)][[1]]}))
}


#--------------------------------------------------------------------------------
#STILL WORKING ON THESE BELOW...
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
orf2entrez <- function(x) {
  d = org.Sc.sgdENTREZID
  mapped_probes <- mappedkeys(d)
  dd = as.list(d[mapped_probes])
  genelist = sapply(x, function(x) {dd[as.character(x)][[1]]})
  unlist(genelist)
  genelist = sapply(genelist, function(x) {
    if (nchar(x) != 7) {
    return(paste(rep("0",7-nchar(x)),x))
  } 
  })
  genelist = paste("GO:",genelist,sep="")
  return(genelist)
}

#--------------------------------------------------------------------------------
goHyperTest <- function(genelist, geneuniverse) {
  frame = toTable(org.Sc.sgdGO)
  goframeData = data.frame(frame$go_id, frame$Evidence, frame$systematic_name)
  goFrame = GOFrame(goframeData, organism = "Saccharomyces cerevisiae")
  goAllFrame = GOAllFrame(goFrame)
  gsc = GeneSetCollection(goAllFrame,setType=GOCollection())
  
  params = GSEAGOHyperGParams(name = "Saccharomyces cerevisiae GO",
                              geneSetCollection = gsc, geneIds = genelist, universeGeneIds = geneuniverse,
                              ontology="BP",pvalueCutoff=0.05,conditional=FALSE,testDirection="over")
  over = hyperGTest(params)
  return(over)
}

TestCor <- function(eig,dat){
  m = apply(dat,1,cor,y=eig,method="pearson")
  return(m)
}

#--------------------------------------------------------------------------------
EigCor <- function(eigengenes,data,ntotest) {
  corr.mat = matrix(0,nrow(data),ntotest)
  for (i in 1:ntotest){
    eig = eigengenes[i,]
    corr.mat[,i] = TestCor(eig,data)
  }
  rownames(corr.mat) = rownames(data)
  return(corr.mat)
}

#--------------------------------------------------------------------------------
findThresh <- function(dat,eigengenes,neigen,thresh,nrep) {
  p.corr=NULL
  for (i in 1:nrep){
    if(i%%10==0){
      print(i)
    }
    p.mat = matrix(sample(dat,size=length(dat),replace=F),nrow(dat),ncol(dat))
    p.corr = c(p.corr,EigCor(eigengenes,p.mat,neigen))
  }
  return(quantile(p.corr,thresh))
}

#--------------------------------------------------------------------------------
plotEigs <- function(eigengene,neigen,eigenmatch) {
  match.eig1 = eigenmatch
  eig.melt = data.frame(eigengene=neigen,eigengenes[neigen,],variable=rep(c(0,5,10,15,20,30,45,60),each=1))
  colnames(eig.melt)[2] = "value"
  
  match.eig1.m = data.frame(orf=rownames(match.eig1), match.eig1)
  match.eig1.m = melt(match.eig1, idvars=orf)
  match.eig1.m$Var2=rep(c(0,5,10,15,20,30,45,60),each=nrow(match.eig1))
  p = ggplot()+geom_line(aes(x=Var2,y=value,group=Var1),match.eig1.m,alpha=0.5,colour="grey")+theme_bw()+
    geom_line(aes(x=variable,y=value),eig.melt,colour="red",size=1)+xlab("time")+ylab("fold change")
  print(p)
}