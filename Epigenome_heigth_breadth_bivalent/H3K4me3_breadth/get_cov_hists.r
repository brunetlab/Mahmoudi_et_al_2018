calibrate_hists <- function(my.prefix,my.path,my.ref.hist,my.ref.name) {
  my.pattern <- paste("hist_",my.prefix,"(.)+bed$",sep="")
  my.files <- list.files(path = my.path, pattern = my.pattern, full.names=TRUE)
  
  my.order <- c(12,2:11,13:20,1)
  
   density.per.peak.list <- vector(mode="list",length=length(my.files))

   for (i in 1:length(my.files)) {
     density.per.peak.list[[i]] <- read.table(my.files[my.order[i]],header=F,sep="\t")
   }
  
  
  my.palette <- rainbow(20)
  my.names <- rep(0,20)
  for (i in 1:20 ){
    index <- seq(5,100,5)[i]
    my.names[i] <- paste(index,"% depth of ", my.prefix,sep="")
  }
  
  my.res.ix <- get_max_ks_P(density.per.peak.list,my.ref.hist);
    
  my.pdf.name <- paste("hists_of_coverage_of_REF_peaks_",my.prefix,".pdf",sep="")
    
  pdf(my.pdf.name)
  plot(0:39,my.ref.hist[1:40],type='l',lwd=2,ylim=c(0,1),ylab="Histogram of H3K4me3 reference peak coverage",xlab="Coverage (X)",main=my.prefix)
  for (i in 1:length(my.files)) {
    lines(0:39,density.per.peak.list[[i]]$V5[1:40],type='l',col=my.palette[i])
  }
  text(5,1,paste(seq(5,100,5)[my.res.ix]," % downsampling",sep=""),cex = 0.5)
  legend("bottomright",c(my.ref.name,my.names),col=c("black",my.palette),bty='n',pch='_',cex = 0.65 ,pt.cex=1.5)
  dev.off()
  
  plot(0:39,my.ref.hist[1:40],type='l',lwd=2,ylim=c(0,1),ylab="Histogram of H3K4me3 reference peak coverage",xlab="Coverage (X)",main=my.prefix)
  for (i in 1:length(my.files)) {
    lines(0:39,density.per.peak.list[[i]]$V5[1:40],type='l',col=my.palette[i])
  }
  text(5,1,paste(seq(5,100,5)[my.res.ix]," % downsampling",sep=""),cex = 0.5)
  legend("bottomright",c(my.ref.name,my.names),col=c("black",my.palette),bty='n',pch='_',cex = 0.65 ,pt.cex=1.5)
  
  print(paste(seq(5,100,5)[my.res.ix]," % downsampling, index : ",my.res.ix,sep=""))
  density.per.peak.list

}

## max p-value gives closest ECDF
get_max_ks_P <- function(density.per.peak.list,my.ref.hist) {
  my.ks.ps <- rep(0,length(density.per.peak.list))
  for (i in 1:length(density.per.peak.list)) {    
    tmp.vec <- ks.test(density.per.peak.list[[i]]$V5, my.ref.hist, alternative = "two.sided")
    my.ks.ps[i] <- tmp.vec$p.value
  }
  my.closest <- which(my.ks.ps == max(my.ks.ps))
  my.closest
}
