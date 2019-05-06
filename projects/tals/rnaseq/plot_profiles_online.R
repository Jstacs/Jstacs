library(scales)
args=c("BXOR1_RVDs/Predicted_binding_sites_for_TalAD9/Profile_for_Chr1_24285654_-.tsv",2,3,"AD9","-","B8-12","24285654","Chr1")
if (length(args)<8) {
  stop("Please submit the following arguments: ProfileFile, numberOfTreatmentReplicates, numberOfControlReplicates, NameOfTALE, Strand, Strain, PositionOfBindingSite", "Chromosome", call.=FALSE)
}else{
  file<-args[1]
  numTreatment<-as.numeric(args[2])
  numControl<-as.numeric(args[3])
  colsTreatment<-seq(2,numTreatment+1,1)
  colsControl<-seq(numTreatment+2,numTreatment+1+numControl,1)
  TALE<-args[4]
  strand<-args[5]
  strain<-args[6]
  BSPos<-as.numeric(args[7])
  chr<-args[8]
  chrNr<-as.numeric(substring(chr,4))
}

  gff<-read.table("all.gff3",sep="\t",quote="",stringsAsFactors = F)
  x<-read.table(file,sep="\t",header=T)
  
  pdf(gsub(x = file,pattern = "\\.tsv",replacement = ".pdf"),width = 15,height = 5*0.85);
  sub<-gff[gff[,1]==chr & gff[,5]>x[1,1] & gff[,4]<x[nrow(x),1] & gff[,3]%in%c("mRNA","exon"),]
  index_mRNA<-vector()
  if(nrow(sub)>0){
    index_mRNA<-which(sub[,3]=="mRNA")
  }  
  
  
  if(length(index_mRNA)>0){
    yoff<-ceiling(max(x[c(colsTreatment,colsControl)]))*0.1*length(index_mRNA)
    index_mRNA<-which(sub[,3]=="mRNA")
    yoff_mRNA<--(seq(1,length(index_mRNA)*1,1)*yoff)-0.3*yoff
    plot(1,axes=FALSE,col=rgb(0.35,0.7,0.9),t="l",ylim=range(x[,c(colsTreatment,colsControl)]) - c(yoff*(length(index_mRNA)+1.5),-0.5),xlim=c(min(x[,1]),max(x[,1])),xlab="genomic position",ylab="",cex.lab=1.3)
  }else{
    plot(1,axes=FALSE,col=rgb(0.35,0.7,0.9),t="l",ylim=range(x[,c(colsTreatment,colsControl)]) - c(diff(range(x[,c(colsTreatment,colsControl)]))*0.1,0),xlim=c(min(x[,1]),max(x[,1])),xlab="genomic position",ylab="")
  }
  
  axis(side=1, cex.axis=1.3)
  axis(side=2,0:ceiling(max(x[c(colsTreatment,colsControl)])), cex.axis=1.3)
  mtext(chr,1,adj=-0.05,line=1,cex=1.3)
  mtext(strain,2,line=2.5,col = "gray60",at=ceiling(max(x[c(colsTreatment,colsControl)]))/2,cex=1.5)
  
  if(nrow(sub)>0){
    xoff<-(apply(cbind(sub[sub[,3]=="mRNA",4],x[1,1]),1,max) + apply(cbind(sub[sub[,3]=="mRNA",5],x[nrow(x),1]),1,min) )/2
    idx_assign<-data.frame(index_mRNA,yoff_mRNA)
    
    if(sum(sub[,3]=="mRNA"&sub[,7]=="+")>0){
      akt_yoff=idx_assign[index_mRNA %in% which(sub[,3]=="mRNA"&sub[,7]=="+"),2]
      ax0<-sub[sub[,3]=="mRNA"&sub[,7]=="+",4]
      ax1<-sub[sub[,3]=="mRNA"&sub[,7]=="+",5]
      arrows(x0 = ax0,y0 = akt_yoff,x1 = ax1,y1 = akt_yoff)
    }
    
    if(sum(sub[,3]=="mRNA"&sub[,7]=="-")>0){
      
      akt_yoff=idx_assign[index_mRNA %in% which(sub[,3]=="mRNA"&sub[,7]=="-"),2]
      ax0<-sub[sub[,3]=="mRNA"&sub[,7]=="-",5]
      ax1<-sub[sub[,3]=="mRNA"&sub[,7]=="-",4]
      arrows(x0 = ax0,y0 = akt_yoff,x1 = ax1,y1 = akt_yoff)
    }
    
    text(x = xoff, y = idx_assign[,2],pos=1,labels = gsub("ID=","",gsub(";.*","",sub[sub[,3]=="mRNA",9])),cex=1.2)
    
    if(sum(sub[,3]=="exon")>0){
      idx_exon<-which(sub[,3]=="exon")
      akt_yoff<- vector()
      m<-1
      for(j in 1:length(sub[,1])){
        if(sub[j,3]=="exon"){
          akt_yoff<-c(akt_yoff,idx_assign[m-1,2])
        }
        else if(sub[j,3]=="mRNA"){
          m<-m+1
        }
      }
      ax0<-sub[sub[,3]=="exon",4]
      ax1<-sub[sub[,3]=="exon",5]
      segments(x0 = ax0,y0 = akt_yoff,x1 = ax1,y1 = akt_yoff,col = 4,lwd = 4,lend=1)
    }
  }
  
  ma<-max(x[,c(colsTreatment,colsControl)])
  lastCol<-max(colsControl)+1
  cols<-ifelse(x[-1,lastCol],4,1)
  streches<-x[which(x[,lastCol]),1]
  if(length(streches)>0){
    starts <- min(streches)
    ends<-vector()
    for(i in 1:(length(streches)-1)){
      if(!(streches[i]==(streches[i+1]-1))){
        ends<-c(ends,streches[i])
        starts<-c(starts,streches[i+1])
      }
    }
    ends<-c(ends,max(streches))
  }
  
  posTALE=min(x[,1])+(max(x[,1])-min(x[,1]))/2
  abline(v=posTALE, col=4, lwd=3)
  mtext(TALE,3,line=-0.8,at=c(posTALE+20), col = "gray60",adj=c(0),cex=1.3)
  if(strand=="-"){
    mtext("<-",3,line=0.2,at=c(posTALE-70), col = 4,adj=c(0),cex=1.3)
  }
  if(strand=="+"){
    mtext("->",3,line=0.2,at=c(posTALE-20), col = 4,adj=c(0),cex=1.3)
  }

  for(i in 1:(length(starts))){
    rect(starts[i],0,ends[i],max(x[,c(colsTreatment,colsControl)]),col=alpha(4, 0.2),density=20)
  }
  
  for(i in 1:length(colsTreatment)){
    lines(x[,1],x[,colsTreatment[i]],col=rgb(0.35,0.7,0.9))
  }
  for(i in 1:length(colsControl)){
    lines(x[,1],x[,colsControl[i]],col=rgb(0.8,0.4,0))
  }
  
  meanTreatment=apply(x[,colsTreatment], 1, mean)
  meanControl=apply(x[,colsControl], 1, mean)
  
  lines(x[,1],meanTreatment,col=rgb(0.35,0.7,0.9),lwd=3)
  lines(x[,1],meanControl,col=rgb(0.8,0.4,0),lwd=3)
  
  dev.off()






