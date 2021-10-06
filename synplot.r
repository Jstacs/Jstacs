getCumSum <- function(a,i,j,target=sort(as.character(unique(a[,i])))) {
	maxT = sapply( target, function(t) max(as.numeric(a[which(a[,i]==t),j])) )
	cst=c(0,cumsum(maxT))
	names(cst)=c(names(maxT),"sum");
	return( cst );
}

orderContigs <- function( a, i=1, j=3, k=2, uni=unique(a[,i]), minPerc=0.1 ) {
	stat=sapply( uni, function(u) {
		idx=which(a[,i]==u);
		if( length(idx)<25 ) return( NULL )
		v=unique(a[idx,j]);
		if( length(v)> 1 ) {
			res=sapply(v,function(vv) {
				h=which(a[idx,j]==vv)
				pos=as.numeric(a[idx[h],k])
				return( c( length(pos), mean(pos) ) );
			});
			res=res[,order(res[2,])]
			perc=res[1,]/length(idx)
			help=which(perc>=minPerc);
			if( length(help)>=1 ) {
				v=colnames(res)[help]
			} else {
				v=NULL;
			}
		}
		return(v);
	} )
	ul=unlist(stat);
	#my unique
	preordered=c();
	for( z in seq_along(ul) ) {
		if( !(ul[z] %in% preordered) ) {
			preordered=c(preordered,ul[z]);
		}
	}
	
	names(preordered)=NULL;
	return( c(preordered,setdiff(unique(a[,j]),preordered)) );
}

addAxis <- function( side, cs ) {
	n=seq_along(cs)[-1];
	axis(side,(cs[n-1]+cs[n])/2/1E6,names(cs)[n-1],las=2, cex.axis=0.5, col=0, col.ticks=1);
}

synPlot <- function( syn, column, k, target, pch=16, cex=0.1, order=T ) {	
	name=colnames(syn)[column];
	idx=which(syn[,column]!="")
	orientContig = apply(syn[idx,c(3,column)],1,function(r) {
		m=matrix(unlist(strsplit(strsplit(r[2],"; ")[[1]],",")),ncol=5,byrow=T);
		paste(
			as.numeric(r[1])*as.numeric(m[,3]),
			m[,1],
			sep=":"
		)
	} );
	ord = sapply(syn[idx,column],function(r) as.numeric(matrix(unlist(strsplit(strsplit(r,"; ")[[1]],",")),ncol=5,byrow=T)[,4]));
	middle = sapply(syn[idx,column],function(r) matrix(unlist(strsplit(strsplit(r,"; ")[[1]],",")),ncol=5,byrow=T)[,2]);
	gene = sapply(syn[idx,column],function(r) matrix(unlist(strsplit(strsplit(r,"; ")[[1]],",")),ncol=5,byrow=T)[,5]);
	
	all=c();
	for( i in 1:(length(idx)-k) ) {
		if( length(unique(syn[idx[i:(i+k-1)],1])) == 1 ) {
			res=c();
			for( ii in i:(i+k-1) ) {
				res=c(res,unique(orientContig[[ii]]));
			}
			t=table(res);
			index=which(t>=k);
			if( length(index)>0 ) {
				for( nam in names(t[index]) ) {
					dir=as.numeric(strsplit(nam,":")[[1]][1])
					col=ifelse(dir==1,1,2)
					hh=which(orientContig[[i]]==nam)
					for( h in hh ) {
						res = c(syn[idx[i],1],syn[idx[i],2],syn[idx[i],5],strsplit(nam,":")[[1]][2],middle[[i]][h],gene[[i]][h],ord[[i]][h],col);
						start=ord[[i]][h];
						ii=1;
						while( ii < k && ((start+ii*dir) %in% ord[[i+ii]][which(orientContig[[i+ii]]==nam)]) ) {
							help=which(orientContig[[i+ii]]==nam);
							z=which( start+ii*dir == ord[[i+ii]][help] )
							res = rbind(res,c(syn[idx[i+ii],1],syn[idx[i+ii],2],syn[idx[i+ii],5],strsplit(nam,":")[[1]][2],middle[[i+ii]][help[z]],gene[[i+ii]][help[z]],ord[[i+ii]][help[z]],col));
							ii=ii+1;
						}
						if( ii == k ) {
							all=rbind(all,res);
						}
					}
				}
			}
		}
	}
	a=unique(all);
	
	csr = getCumSum( do.call( rbind, sapply(syn[idx,column],function(r) matrix(unlist(strsplit(strsplit(r,"; ")[[1]],",")),ncol=5,byrow=T)[,c(1,2)], USE.NAMES = F) ), 1, 2 )
	if( !order ) {
		cst = getCumSum(syn,1,2)
		#csr = getCumSum( do.call(rbind,sapply(syn[idx,column],function(r) matrix(unlist(strsplit(strsplit(r,"; ")[[1]],",")),ncol=5,byrow=T)[,c(1,2)]) ), 1, 2 )
	} else {
		ord = orderContigs(a,i=4,j=1,k=5, uni=names(csr)[-length(csr)]);
		cst = getCumSum(syn,1,2, target=c(ord,setdiff(unique(syn[,1]),ord)));
		#csr = getCumSum( do.call(rbind,sapply(syn[idx,column],function(r) matrix(unlist(strsplit(strsplit(r,"; ")[[1]],",")),ncol=5,byrow=T)[,c(1,2)]) ), 1, 2, target=orderContigs(a) )
	}
	
	ct=as.vector(a[,1])
	t=as.numeric(a[,2])
	cr=as.vector(a[,4])
	r=as.numeric(a[,5])
		
	x=(cst[ct]+t)/1E6;
	y=(csr[cr]+r)/1E6;
	
	plot(x,y,col=as.numeric(a[,8]),pch=pch,cex=cex,ylim=c(0,csr["sum"])/1E6,xlim=c(0,cst["sum"])/1E6,ylab="",xlab="",axes=F, main=paste(target, " vs ", name, " (k=",k,")",sep=""))
	abline(v=cst/1E6,lty=2,col=gray(0.8),lwd=0.1)
	abline(h=csr/1E6,lty=2,col=gray(0.8),lwd=0.1)
	addAxis(1,cst)
	addAxis(2,csr)
	
	return(a);
}

###test parameter
#inputName="./reference_gene_table.tabular";
#k=5;
#outputName="GSD";

#user parameters
args <- commandArgs(trailingOnly = TRUE)
inputName=args[1]
k=as.numeric(args[2]);
outputName=args[3];

#read table
syn=as.matrix(read.delim(inputName,sep="\t",header=T));

#plot and write syntheny
pdf(paste(outputName,"-synteny-k=",k,".pdf",sep=""))
res=sapply( 7:ncol(syn), function(j) {
	res=synPlot( syn, j, k, name );
	write.table(res,paste(name,"-synteny-",j,"-k=",k,".tabular",sep=""),col.names=F,row.names=F,sep="\t", quote=F)
} )
dev.off();
