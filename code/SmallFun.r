getIBS<-function(geno,weights=1)
{
	n<-nrow(geno);
	m<-ncol(geno);
	gtemp1<-geno;
	gtemp2<-geno;
	gtemp1[geno==2]<-1;
	gtemp2[geno==1]<-0;
	gtemp2[geno==2]<-1;
	gtemp<-cbind(gtemp1,gtemp2);
	Inner<-gtemp%*%diag(weights,nrow=2*m,ncol=2*m)%*%t(gtemp);
	X2<-matrix(diag(Inner),nrow=n,ncol=n);
	Y2<-matrix(diag(Inner),nrow=n,ncol=n,byrow=T);
	Bound<-sum(matrix(weights,nrow=1,ncol=m)*2);
	IBS<-Bound-(X2+Y2-2*Inner);
	I1=diag(nrow(IBS));
	J1=matrix(1/nrow(IBS),nrow=nrow(IBS),ncol=ncol(IBS));
	IBSnew=(I1-J1) %*% IBS %*% (I1-J1);
	diag(IBSnew)=0
	IBSnew;
}

sampleChr<-function(input,nloc)
{
	load(input);
	Chr=CHR
	pos=as.numeric(as.matrix(Chr$pos[,2]));
	geno=Chr$geno;
	Genes=list();
	for(i in 1:nloc)
	{
		aa=sample(pos,1)
		upper=aa+10000
		lower=aa-10000		
		select=which(pos<upper & pos>lower);
		Genes[[i]]=colnames(geno)[select]
		geno=geno[,-select];
		pos=pos[-select]
	}
	rm(CHR);gc();
	result=list();
	result$Gene=Genes;
	result$geno=Chr$geno[,colnames(Chr$geno) %in% unlist(Genes)];
	result;
}

sampleGeno<-function(chrs,nloc,ndrop=0)
{
	RESULT=list();
	GenoDrop=NULL;
	for(i in 1:length(chrs))
	{
		RESULT[[i]]=sampleChr(chrs[i],nloc[i])
	}
	
	Genes=list();
	Geno=NULL;
	num=1;
	for(i in 1:length(chrs))
	{
		if(num==1)Geno=RESULT[[i]]$geno
		if(num>1) Geno=cbind(Geno,RESULT[[i]]$geno)
		for(j in 1:nloc[i])
		{
			if(ndrop[i]==0)
			{
				Genes[[num]]=RESULT[[i]]$Gene[[j]]
				num=num+1;
			}
			if(ndrop[i]!=0)
			{
				if(ndrop[i]>1)
				{
					print("More than 1 gene is dropped, switch to 1");
					ndrop[i]=1
				}
				if(j!=nloc[i])
				{
					Genes[[num]]=RESULT[[i]]$Gene[[j]]
					num=num+1;
				}
				if(j==nloc[i])GenoDrop=RESULT[[i]]$Gene[[j]]
			}		
		
		}
	}

		
	result=list();
	result$Genes=Genes;
	result$Geno=Geno[,!colnames(Geno) %in% GenoDrop];
	result$Gdrop=Geno[,colnames(Geno) %in% GenoDrop];
	result;
}


Minor.Afreq<-function(geno)
{
	if(is.null(dim(geno)))
	{
		Afreq=mean(geno)/2
	}
	if(!is.null(dim(geno)))
	{
		Afreq=colMeans(geno,na.rm = T)/2
	}
	t1=Afreq>0.5;
	Afreq[t1]=1-Afreq[t1];
	Afreq;
}
uniform.fmatrix<-function(X) #family distance matrix, with diag = 0
{
  n<-length(X[,1])
  temp<-matrix(0,n,n)
  #checking parent-kid
  temp<-temp+outer(X[,1],X[,1],Check.share)
  #diag(temp)<-0
  return(temp)
}
Check.share<-function(x,y){return(as.numeric(x==y))}



