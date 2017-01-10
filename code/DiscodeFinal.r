## Modify the error parts, as if the family is not sorted the codes are not right ##
ERROR.ENV<-function(wtinput,Pop)
{
	error=rep(0,nrow(Pop))
	for (j in unique(Pop[,1]))
	{
		t1=Pop[,1]==j
		error[t1]=getErrorEnv(Pop[Pop[,1]==j,],wt=wtinput)
	}
	error
}


getErrorEnv<-function(Pop,wt,seed=NULL)
{
	Uni<-uniform.fmatrix(Pop[,1:4])#uniform correlation
	Sigma<-Uni*wt #covariance matrix of random effect
	k=ncol(Uni)
	if(length(seed)>0){set.seed(seed);}
	r.effect<-as.matrix(mvrnorm(1,rep(0,k),Sigma))#random effect
	y=r.effect;
	y;
}


DiscodeFinal<-function(chrsinput,nlocsave,nloctrue,effect,ndrop=0,drop.effect=0,lprop=3,sd.err=1,wtinput=1,
pedscheme=2,Noff1=1,Noff2=0,pfamily=2/3,binary=FALSE,scheme="Equal")
{
	samples=sampleGeno(chrsinput,nloc=nlocsave,ndrop=ndrop);
	gdrop=samples$Gdrop;
	genes=samples$Genes;
	geno=samples$Geno;
	geno=as.matrix(geno)
	MAF=Minor.Afreq(geno);
	
	exclude=MAF==0;	
	geno=geno[,!exclude]
	MAF=MAF[!exclude]
	bb=sample(1:nrow(geno),nrow(geno));
	geno=geno[bb,]

	################################
	##### Pedigree generation ######
	################################
	ggtmp=as.matrix(cbind(geno,gdrop))
	pedi=pedigree_gen(pedscheme,Noff1,Noff2,ggtmp);
	GENOSAVE=pedi$GENOSAVE[,1:ncol(geno)];
	PEDSAVE=pedi$PEDSAVE;
	GDROP=pedi$GENOSAVE[,(1+ncol(geno)):ncol(ggtmp)]
	rm(ggtmp)

	################################
	##### Simulation disease #######
	################################
	Y.drop=rep(0,nrow(GENOSAVE));
	if(ndrop!=0) Y.drop=drop.effect*apply(GDROP,1,sum);

	y=rep(0,nrow(GENOSAVE));
	geno=GENOSAVE;
	for(i in 1:nloctrue)
	{
		includetemp=colnames(GENOSAVE) %in% genes[[i]];
		temp=(geno[,colnames(GENOSAVE) %in% genes[[i]]])
		temp=matrix(as.numeric(temp),ncol=ncol(temp))
		maf=MAF[includetemp];
		aa=sample(1:ncol(temp),ncol(temp)/lprop);
		aa=aa[order(aa)];
#		ttgene[[i]]=colnames(GENOSAVE)[colnames(GENOSAVE) %in% genes[[i]]][aa]
		if(scheme=="Equal") wt=rep(1,length(aa))*effect[i]
		if(scheme=="Beta")  wt=dbeta(maf[aa],1,25)^2*effect[i]
		if(scheme=="Log")   wt= -log(maf[aa],base=10)*effect[i]
		if(scheme=="Sqrt")  wt= 1/(maf[aa]*(1-maf[aa]))*effect[i]
		wt=matrix(wt,nrow=nrow(temp),ncol=length(aa),byrow=T);
		y=y+rowSums(temp[,aa]*wt)
	}
	Y.gene=y;

	family=cbind(PEDSAVE[,'FAMID'],PEDSAVE[,'ID'],PEDSAVE[,'DADID'],PEDSAVE[,'MOMID']);
	colnames(family)=c('FAMID','ID','DADID','MOMID');
	error.env=ERROR.ENV(wtinput,Pop=family)

	Y=Y.gene+error.env+Y.drop;

	if(binary)
	{
		Y=Y-mean(Y);
		Y[Y>400]=400;
		Y[Y<(-400)]=-400;
		y=exp(Y)/(1+exp(Y))
		dis=as.numeric(runif(nrow(GENOSAVE),0,1)<y);
	}
	if(!binary)
	{
		dis=Y+rnorm(nrow(GENOSAVE),sd=sd.err)
	}

	result=list();
	result$GENOSAVE=GENOSAVE;
	result$PEDSAVE=PEDSAVE;
	result$DIS=dis
	result$GENE=genes
	result;

}
