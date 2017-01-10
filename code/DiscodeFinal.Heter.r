
DiscodeFinal.Heter<-function(chrsinput,nlocsave,nloctrue,effect,peddir=NULL,ndrop=0,drop.effect=0,lprop=3,sd.err=1,wtinput=1,
pedscheme=2,Noff1=1,Noff2=0,pfamily=2/3,binary=FALSE,scheme="Equal",heter.prop=1,combineGene=TRUE)
{
	samples=sampleGeno(chrsinput,nloc=nlocsave,ndrop=ndrop);
	gdrop=samples$Gdrop;
	if(combineGene)
	{
		tmp=NULL;
		for(i in 1:length(samples$Genes))
		{
			tmp=c(tmp,samples$Genes[[i]])
		}
		samples$Genes=NULL;
		samples$Genes[[1]]=tmp
	}
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
	pedi=pedigree_gen(pedscheme,Noff1,Noff2,ggtmp,peddir);
	GENOSAVE=pedi$GENOSAVE[,1:ncol(geno)];
	PEDSAVE=pedi$PEDSAVE;
	GDROP=pedi$GENOSAVE[,(1+ncol(geno)):ncol(ggtmp)]
	rm(ggtmp)

	################################
	##### Simulation disease #######
	################################
	Y.drop=rep(0,nrow(GENOSAVE));
	if(ndrop!=0) Y.drop=drop.effect*apply(GDROP,1,sum);


	family=cbind(PEDSAVE[,'FAMID'],PEDSAVE[,'ID'],PEDSAVE[,'DADID'],PEDSAVE[,'MOMID']);
	colnames(family)=c('FAMID','ID','DADID','MOMID');
	error.env=ERROR.ENV(wtinput,Pop=family)

	y=rep(0,nrow(GENOSAVE));
	geno=GENOSAVE;
	for(i in 1:nloctrue)
	{
		includetemp=colnames(GENOSAVE) %in% genes[[i]];
#################################################################################
########### Adding heterogeneity ##############		
##################################################################################
		temp=(geno[,colnames(GENOSAVE) %in% genes[[i]]])
		temp=matrix(as.numeric(temp),ncol=ncol(temp))
		maf=MAF[includetemp];
		aa=sample(1:ncol(temp),ncol(temp)/lprop);
		aa=aa[order(aa)];
##		ttgene[[i]]=colnames(GENOSAVE)[colnames(GENOSAVE) %in% genes[[i]]][aa]
		if(scheme=="Equal") wt=rep(1,length(aa))*effect[i]
		if(scheme=="Beta")  wt=dbeta(maf[aa],1,25)^2*effect[i]
		if(scheme=="Log")   wt= -log(maf[aa],base=10)*effect[i]
		if(scheme=="Sqrt")  wt= 1/(maf[aa]*(1-maf[aa]))*effect[i]
		wtsave=wt
		wt=matrix(wt,nrow=nrow(temp),ncol=length(aa),byrow=T);
		if(heter.prop==1) y=y+rowSums(temp[,aa]*wt)

		n.family=unique(family[,'FAMID']);

		if(heter.prop!=1)
		{
			wtt=matrix(0,nrow=nrow(temp),ncol=length(aa))
			if(heter.prop==0) n.group=length(n.family);
			if(heter.prop!=0) n.group=round(1/heter.prop);
			if(n.group > length(aa)) n.group=length(aa);
		#	n.loci.each=floor(length(aa)/n.group);
			geno.num.seq=floor(seq(from=1,to=length(aa),length.out=(n.group+1)))
			fam.num.seq=floor(seq(from=1,to=length(n.family),length.out=(n.group+1)))

			for(j in 1:n.group)
			{
				fam=n.family[fam.num.seq[j]:(fam.num.seq[j+1]-1)]
				wt.tmp=rep(0,length(aa))
				wt.tmp[geno.num.seq[j]:(geno.num.seq[j+1]-1)]=wtsave[geno.num.seq[j]:(geno.num.seq[j+1]-1)]
				wtt[family[,'FAMID'] %in% fam,]=matrix(wt.tmp,nrow=sum(family[,'FAMID'] %in% fam),ncol=length(aa),byrow=T);
			}
	
			y=y+diag(temp[,aa] %*% t(wtt))
		}
	}
	Y.gene=y;


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
