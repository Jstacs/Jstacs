#This R scricpt allows to visualize the results of the module Analyzer if w=YES

catCol=c("yellow","blue")
#functions
plotCategorie <- function( bf1, main="true vs. predicted transcripts (best)" ) {
	categories = c(
		length(which(bf1==0)),
		length(which(bf1>0 & bf1<1)),
		length(which(bf1==1)),
		length(which(bf1<0))
	)
	
	#pie( categories, labels=c("no overlap", "any overlap", "identical", "alternative"), col=c("red","yellow","green","orange"), main=main )
	
	barplot( categories, names.arg=c("no overlap", "any overlap", "identical", "alternative"), col=c("red","yellow","green","orange"), main=main, las=2, cex.names=0.8 )
	p=seq(0,1,by=0.1)
	all=sum(categories);
	axis(4,at=all*p,p*100,col="gray",las=1)
	abline(h=seq(0,1,by=0.05)*all,lty=2,col="gray",col.axis="gray");
	
	return( categories/sum(categories) );
}

getRelevant <- function( x ) {
	u=unique(x)
	selected=unlist(sapply(u,function(h){ idx=which( x==h ); if( length(idx)==1 ) { return(idx); } else { return(idx[c(1,length(idx))]); } }))
	return( selected );
}

combiPlot <- function( t, names, all=sum(t), col=catCol ) {
	max=max(t)
	y=seq(-max,-0.1*max,length=length(names));
	
	plot(seq_along(t),t,type="h",axes=F,ylab="",xlab="",xlim=c(0.5,length(t)+0.5),ylim=max*c(-1,1),lwd=2);
	axis(2,seq(0,max,length=5),las=1)
	axis(4,at=seq(0,1,by=0.05)*all,seq(0,100,by=5),las=1,col="gray",col.axis="gray");
	abline(h=seq(0,1,by=0.05)*all,lty=2,col="gray");
	
	h=sapply( names(t), function(n) as.integer(strsplit(n,"")[[1]]) );
	image(t(h),x=seq_along(t),y=y,add=T,axes=F,col=col);
	axis(2,y,gsub("\\."," ",names),cex=0.5,las=1)
}

f1VsSn <- function( fName, add=F, col=1 ) {
	data=read.delim(fName,sep="\t")

	bf1=data[,"transcript.bestF1"]
	f1=data[,"transcript.F1"]
	pIdx=which(!is.na(bf1) & (bf1<=0 | bf1==f1))
	
	numPred=length(which(data[,"predictionId"]!=""))
	
	##all
	data2=data[pIdx,]

	bf1=data2[,"transcript.bestF1"]
	f1=data2[,"transcript.F1"]
	presel=f1;
	x = rev(seq_along(presel))/length(presel)*100;
	
	idx=which(is.na(presel))
	
	g=presel;
	g[idx]=0;
	sf1_1 = sort(g)
	selected1=getRelevant( sf1_1 );
	
	h=presel;
	h[idx]=abs(bf1[idx])
	sf1_2 = sort(h)
	selected2=getRelevant( sf1_2 );
	
	if( !add ) {
		plot(x[selected1],sf1_1[selected1],ylab="best transcript F1",xlab="transcript sensitivity",type="l", main="true transcript vs. best prediction", col=col)
		lines(x[selected2],sf1_2[selected2],col=col,lty=2)
	} else {
		lines(x[selected1],sf1_1[selected1],col=col);
		lines(x[selected2],sf1_2[selected2],col=col,lty=2)
	}
	abline(v=c(0,100),col=gray(0.8),lty=2)
	abline(h=0:1,col=gray(0.8),lty=2)
	
	perfect = length(which(g==1))
	sn = perfect / length(g) *100;
	prec = perfect / numPred *100
	
	return( c(sn, prec) );	
}

analyze <- function( data, info, numPred,col=catCol ) {
	initial=par();
	plot.new()
	text(0.5,0.5,info,cex=2,col=2)

	#best f1 categories
	bf1=data[,"transcript.bestF1"]
	plotCategorie(bf1);
	
	#plot best transcript F1 vs. sensitivity
	bf1=data[,"transcript.bestF1"]
	f1=data[,"transcript.F1"]
	
	presel=f1;
	x = rev(seq_along(presel))/length(presel)*100;
	
	idx=which(is.na(presel))
	
	g=presel;
	g[idx]=0;
	sf1 = sort(g)
	selected=getRelevant( sf1 );
	plot(x[selected],sf1[selected],ylab="best transcript F1",xlab="transcript sensitivity",type="l", main="true transcript vs. best prediction")
	
	h=presel;
	h[idx]=abs(bf1[idx])
	sf1 = sort(h)
	selected=getRelevant( sf1 );
	lines(x[selected],sf1[selected],col=4)
	
	perfect = length(which(g==1))
	sn = perfect / length(g) *100;
	prec = perfect / numPred *100
	abline(v=sn,lty=3,col=2,lwd=2)
	text(sn,1,paste("sensitivity = ",round(sn,1),"\nprecision = ",round(prec,1),sep=""),col=2,adj=c(1.1,1.2))
	
	abline(v=100-length(which(h==0)) / length(g) *100,lty=3,col=3,lwd=2)
	
	abline(v=c(0,100),col=gray(0.8),lty=2)
	abline(h=0:1,col=gray(0.8),lty=2)
	legend("bottomleft", fill=c(1,4), c("conservative (i.e. bestF1>=0)","liberal (i.e. abs(bestF1))"), bg=gray(0.999) );

	
	#analyze any overlap
	data3=data[which(f1>0 & f1<1),]
	##TODO further subselection
	
	##single
	plot.new()
	text(0.5,0.5,"any overlap",cex=2)
	text(0.5,0.3,"single stat")
	
	#histogram/ecdf
	name="feature.difference";
	plot( ecdf(data3[,name]),main=gsub("\\."," ",name))
	
	#ignore "tp", "fn", "fp", "perfect.features"
	
	#categories true, false
	namesBoolean=c("start","end","first.feature","last.feature")
	#categories: 0, >0
	namesZero=c("additional.upstream.features.truth", "additional.internal.features.truth", "additional.downstream.features.truth", "additional.upstream.features.prediction", "additional.internal.features.prediction", "additional.downstream.features.prediction", "intron.retention.in.truth", "intron.retention.in.prediction")
	namesOne=c("donor","acceptor");
	names=c(namesBoolean,namesZero,namesOne);
	val=c(rep("true",length(namesBoolean)),rep("0",length(namesZero)),rep("1",length(namesOne)));
	#alt=c(rep("false",length(namesBoolean)),rep("at least 1",length(namesZero)));
	
	oldPar=par();
	mat = sapply(seq_along(names),function(i) ifelse(!is.na(data3[,names[i]]) & data3[,names[i]]==val[i],0,1) )
	v=apply(mat,2,function(x) { s=sum(x)/length(x); return( c(s,1-s) ) })
	par(mar=c(18,4,1,1)+0.1)
	o=order(v[1,])
	barplot(100*v[,o],col=rev(col),names.arg=gsub("\\."," ",names[o]),las=2)
	par(oldPar);
	
	##complex
	plot.new()
	text(0.5,0.5,"any overlap",cex=2)
	text(0.5,0.3,"combined stat")
	
	par(mar=c(0,18,0,2)+0.1)
	t=table(apply(mat,1,paste,collapse=""))
	#combiPlot(t,names)
	combiPlot(sort(t),names)
	
	mostRelevant=which(t>=0.05*sum(t))
	tt=sort(t[mostRelevant])
	combiPlot(tt,names,sum(t))
	par(initial)
	
	return( c(length(g),numPred,perfect) );
}

myVioplot <- function( res, main, col ) {
	min=min(unlist(lapply(res,min,na.rm=T)));
	max=max(unlist(lapply(res,max,na.rm=T)));
	if( max-min>0 ) {
		#boxplot(res,col=col,main=main,log=ifelse(min<=0,"","y"))
	
		l=length(res)
		plot(0,0,col=0,xlim=c(1,l*2+1),ylim=c(min,max),main=main,log=ifelse(min<=0,"","y"),axes=F,ylab="",xlab="")
		axis(2);
		for( i in 1:l ) {
			idx=which(!is.na(res[[i]]))
			if( length(idx)>0 ) {
				u=unique(res[[i]][idx]);
				if( length(u)==1 ) {
					#TODO
					points(i*2,u,pch="+",col=col[i]);
				} else {
					vioplot(res[[i]][idx],at=i*2,col=col[i],add=T)
				}			
			}
		}
		axis(1,2*(1:l),names(res));
		return(T);	
	}
	return(F)
}

minus = function(x,y) { return( x-y ) }
div = function(x,y) { return( x/y ) }
	

#TODO user input from command line
input=list();
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	stop("Missing argument.", call.=FALSE)
	
	#old test case
	#input[[1]]="E:/Analyzer/GeMoMa/4ref-neu/neu/Comparison_0_1.tabular"
	#
	#input=c("4ref_all=E:/Analyzer/GeMoMa/4ref-neu/neu/naive1_all/Comparison_0.tabular",
	#	"4ref_complete=E:/Analyzer/GeMoMa/4ref-neu/neu/naive2_all/Comparison_0.tabular",
	#	"4ref_f=E:/Analyzer/GeMoMa/4ref-neu/neu/naive3_all/Comparison_0.tabular",
	#	"4ref_f,atf=E:/Analyzer/GeMoMa/4ref-neu/neu/naive4_all/Comparison_0.tabular",
	#	"4ref_default=E:/Analyzer/GeMoMa/4ref-neu/neu/Comparison_0_1.tabular")
} else {
	for( i in seq_along(args) ) {
	  	input[[i]]=args[i];
	}
}

if( length(input)==1 ) {
	##analyze one annotation and compare it with thetruth

	#read data
	data=read.delim(input[[1]],sep="\t")
	
	#plots
	pdf(paste(input,"-comparison-with-truth.pdf",sep=""))
	
	##preselection
	bf1=data[,"transcript.bestF1"]
	f1=data[,"transcript.F1"]
	pIdx=which(!is.na(bf1) & (bf1<=0 | bf1==f1))
	
	numPred=length(which(data[,"predictionId"]!=""))
	
	##all
	data2=data[pIdx,]
	analyze(data2,"all annotations",numPred);
	
	##expressed (if data available)
	expressed=c("transcript.tie","transcript.tpc")
	if( all(expressed %in% colnames(data2)) ) {
		tie=data2[,expressed[1]];
		tpc=data2[,expressed[2]];
		exp=which(tpc==1 & (is.na(tie) | tie==1));
		analyze(data2[exp,],"all experimentally verified annotations",numPred);
	}
	dev.off();	
	
	library(vioplot)
	pdf(paste(input[[1]],"-visulize-attributes.pdf",sep=""))
	
	cn=colnames(data);
	idx=which(cn=="prediction.end.position")
	cn=cn[-(1:idx)]
	
	perfectPrediction=which(data[,"transcript.F1"]==1)
	bestOverlapping=which(bf1==f1 & f1<1)
	otherOverlapping=which(bf1!=f1 & f1<1)
	onlyPrediction=which(data[,"transcriptId"]=="")
	col=c(3,4,2,5)
	
	all=list()
	barplot( c(length(perfectPrediction), length(bestOverlapping), length(otherOverlapping), length(onlyPrediction)), col=col, names.arg=c("perfect","best overlapping","other overlapping","only prediction") )
	for( n in cn ) {
		res=list();
		res[["perfect"]]=data[perfectPrediction,n]
		res[["best overlapping"]]=data[bestOverlapping,n]
		res[["other overlapping"]]=data[otherOverlapping,n]
		res[["only prediction"]]=data[onlyPrediction,n]
		h=tryCatch( myVioplot(res, n, col) , error = function(w) { print(paste("Could not plot column",n)); return(F) } )
		if( h ) all[[n]]=res
	}
	
	combi = matrix(NA,ncol=2,nrow=8);
	combi[1,]=c("prediction.aa","prediction.raa")
	combi[2,]=c("prediction.ce","prediction.rce")
	combi[3,]=c("prediction.pAA","prediction.iAA")
	combi[4,]=c("prediction.score","prediction.aa")
	combi[5,]=c("prediction.score","prediction.maxScore")
	combi[6,]=c("prediction.lpm","prediction.raa")
	combi[7,]=c("prediction.maxGap","prediction.raa")
	combi[8,]=c("prediction.score","prediction.bestScore")
	
	
	func = c(minus,minus,minus,div,div,div,div,div);
	op= c("-","-","-","/","/","/","/","/")
	for( n in 1:nrow(combi) ) {
		res=list();
		f=func[[n]];
		res[["perfect"]]=f(data[perfectPrediction,combi[n,1]],data[perfectPrediction,combi[n,2]])
		res[["best overlapping"]]=f(data[bestOverlapping,combi[n,1]],data[bestOverlapping,combi[n,2]])
		res[["other overlapping"]]=f(data[otherOverlapping,combi[n,1]],data[otherOverlapping,combi[n,2]])
		res[["only prediction"]]=f(data[onlyPrediction,combi[n,1]],data[onlyPrediction,combi[n,2]])
		name=gsub("prediction.","",paste(combi[n,1],op[n],combi[n,2]))
		h=tryCatch( myVioplot(res, name, col), error = function(w) { print(paste("Could not plot:",name)); return(F) } )
		if( h ) all[[name]]=res;
	}
	dev.off();
} else {
	#compare different annotations
	pdf("comparison.pdf")
	res = t(sapply(seq_along(input),function(n){
		name=input[[n]];
		split=strsplit(name,"=")[[1]]
		name=split[1];
		fName=split[ifelse( length(split)==1, 1, 2 )];
		return( c(name, f1VsSn(fName,n>1,n)) );
	}))
	colnames(res)=c("name","sn","sp");
	legend("bottomleft", fill=seq_along(input), res[,1], bg=gray(0.999) );

	par(mar=c(5,4,2,4)+0.1)
	m=matrix(as.numeric(t(res[,-1])),nrow=2);
	barplot(m,names.arg=res[,1],beside=T,las=2,ylim=c(0,100),col=catCol)
	F1=apply(m,2,function(x) 2*x[1]*x[2] / (x[1]+x[2]))
	lines(2+3*(seq_along(input)-1), F1, col=2,lty=2,type="b",pch=16)
	a=axTicks(2);
	axis(4,a,seq(0,1,length=length(a)),col=2,col.axis=2,las=2)
	mtext("F1",4,line=3,col=2)
	legend("top", fill=catCol, c("transcript sensitivity", "transcript precision"),ncol=2);
	dev.off()
	
	cbind(res,F1)
}