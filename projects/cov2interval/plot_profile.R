cat("command line arguments: <output-name-PDF> <centromer-file> <reference-name> <coverage_1> ...\n\n")
#user-data
args = commandArgs(trailingOnly=TRUE);
output=args[1];
centromerFile=args[2]
r=args[3]
files=rep(NA,length(args)-3);
for( i in 4:length(args) ) {
	files[i-3]=args[i];
}

cat("output: ",output,"\n",sep="")
cat("centromer info: ",centromerFile," ",r,"\n",sep="")
cat("files: ", paste(files,collapse=", "),"\n",sep="")


#fixed parameters
size=1E6
max=1E9;
xlim=c(0,max)/1E6
pseudo=1e-2;

##functions
readData <- function( fName ) {
	d=as.matrix(read.delim(fName,sep=""),sep="\t",header=F)
	startAll=as.numeric(d[,2])
	endAll=as.numeric(d[,3])
	covAll=as.numeric(d[,4])	
	return( list(chrAll=d[,1], startAll=startAll, endAll=endAll, covAll=covAll) );
}

computeCov = function(chr, chrs, startAll, endAll, covAll, size=5E5) {
	idx=which(chrs==chr);
	if( length(idx)==0 ) {
		return( matrix(NA,ncol=2,nrow=0) );
	}
	start=startAll[idx];
	end=endAll[idx];
	cov=covAll[idx]

	minPos=0;#min(start)
	maxPos=max(end)
	grid=c(seq(minPos,maxPos,by=size)[-1],maxPos+1)
	
	h=which(cov==1);
	idx1=(start[h]-minPos)%/%size + 1
	idx2=(end[h]-minPos)%/%size + 1
	
	res=rep(0,length(grid));
	for( i in seq_along(h) ) {
		if( idx1[i]==idx2[i] ) {
			res[idx1[i]] = res[idx1[i]] + (end[h[i]]-start[h[i]]+1)
		} else {
			#problem idx1[i]+1>idx2[i]
			z=idx1[i]:(idx2[i]-1);
			s=start[h[i]];
			for( a in z ) {
				res[a] = res[a] + (grid[a]-s)
				s=grid[a]
			}
			res[idx2[i]] = res[idx2[i]] + (end[h[i]]-grid[idx2[i]-1]+1)
		}
	}
	
	num=grid-c(minPos,grid[-length(grid)])
	c=res/num * 100;
	return(cbind(grid,c));
}

plotCov <- function( chr, res, col=1, add=F, median=T, main=chr, xlim=NULL, logarithm=F, pseudo=0, yAxesDisp=seq(0,100,by=10) ) {
	x=res[,1]/1E6
	
	y=res[,2];
	
	yAxes=yAxesDisp;
	if( logarithm ) {
		y=log(pseudo+y);
		yAxes=log(pseudo+yAxes);
	}
	ylim=yAxes[c(1,length(yAxes))]

	if( add ) {
		lines(x,y,col=col);
	} else {
		if( is.null(xlim) ) {
			plot(x,y,col=col,type="l",xlab="position in Mb",ylab="coverage>0 in %",main=main, cex.lab=1.25, cex.lab=1.5, cex.main=1.75, ylim=ylim, axes=F);
		} else {
			plot(x,y,col=col,type="l",xlab="position in Mb",ylab="coverage>0 in %",main=main, cex.lab=1.25, cex.lab=1.5, cex.main=1.75, ylim=ylim, xlim=xlim, axes=F);
		}
		axis(1);
		axis(2,at=yAxes,yAxesDisp);
		box();
	}
	if( median && nrow(res)>0 ) {
		m=median(res[,2]);
		mm=m+ifelse(logarithm,1,0);
		abline(h=mm,col=3,lty=3)
		axis(2,mm,round(m,0),col=3,col.axis=3)
	}
}

#read centromers
readInput <- function( file ) {
	len=as.matrix(read.delim(file,sep="\t",header=T))
	names=len[,1];
	a=t(apply(len,1,function(x) as.numeric(x[-1])));
	rownames(a)=names
	colnames(a)=colnames(len)[-1]
	return(a);
}
if( file.exists(centromerFile) ) {
	cen=readInput(centromerFile);
} else {
	cen=NULL;
}

#read data
data=list();
chrs=c();
for( x in files ) {
	data[[x]]=readData( x );
	chrs=unique(c(chrs,data[[x]]$chrAll))
}
chrs=sort(chrs);

#compute per chr
pdf(output,width=21)
for( chr in chrs ) {
	#parse data
	res=list()
	for( x in files ) {
		res[[x]]=computeCov(chr, data[[x]]$chrAll, data[[x]]$startAll, data[[x]]$endAll, data[[x]]$covAll, size=size);
	}

	#create plots	
	for( i in seq_along(files) ) {
		plotCov( chr, res[[files[i]]], col=i, add=i>1, median=F,xlim=xlim, logarithm=T, pseudo=pseudo, yAxesDisp=c(0,0.01,0.1,1,10,100) )
	}
	if( !is.null(cen) ) {
		abline(v=cen[r,chr],col="gray",lty=2)
		text(cen[r,chr],log(pseudo+100),"centromere",srt=90,adj=c(1,0),col="gray",cex=1.5)
	}
}
dev.off()
