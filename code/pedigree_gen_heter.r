pedigree_gen<-function(pedscheme,Noff1,Noff2,geno,peddir)
{
### Considering Subpopulaton information 
### Pedscheme=2, 2 generation case, trio cases, 2 parents, Noff1 offsrpings ####
### Pedscheme=3, 3 genreation case, 2 grandparents, Noff1 offspring for the second generation, Noff2 offsprings for the third generation per second generation ####
	ped=read.csv(peddir,header=T,sep='\t')
	founderindex=ped[,"Paternal.ID"]==0 & ped[,"Maternal.ID"]==0
	geno=geno[founderindex,];
	ped=ped[founderindex,];
	if(pedscheme==2)
	{
		###################################
		#### Need 2 founder per family ####
		###################################
		genoM=genoF=Popindex=NULL;
		uniPop=unique(ped[,'Population'])
		for(j in 1:length(uniPop))
		{
			tmp=geno[ped[,'Population']==uniPop[j],];
			if(nrow(tmp)%%2==1) tmp=tmp[-1,];
			fid=sample(1:nrow(tmp),nrow(tmp)/2)
			genoF=rbind(genoF,tmp[fid,])
			genoM=rbind(genoM,tmp[-fid,])
			Popindex=rbind(Popindex,matrix(rep(uniPop[j],length(fid),ncol=1)))
		}

		geno=rbind(genoF,genoM)
		genoO=NULL;
		FatherID=1:nrow(genoF);
		MotherID=(nrow(genoF)+1):(2*nrow(genoF));
		
 		for(i in 1:Noff1)
		{
			genoO=rbind(genoO,GeneDropCall(genoF,genoM));
		}
		GENOSAVE=rbind(genoF,genoM,genoO);
		IDSAVE=1:nrow(GENOSAVE);
		MOTHERIDSAVE=c(rep(0,nrow(geno)),rep(MotherID,Noff1));
		FATHERIDSAVE=c(rep(0,nrow(geno)),rep(FatherID,Noff1));
		FAMILYIDSAVE=rep(1:(nrow(geno)/2),(2+Noff1));
		SEXSAVE=c(rep(1,nrow(genoF)),rep(2,nrow(genoM)),sample(c(1,2),nrow(genoM)*Noff1,replace=T))
		PEDSAVE=cbind(IDSAVE,FATHERIDSAVE,MOTHERIDSAVE,SEXSAVE,FAMILYIDSAVE);
		colnames(PEDSAVE)=c("ID","DADID","MOMID","SEX","FAMID");
	}


	if(pedscheme==3)
	{
		Nind=2+Noff1;

		
		genoM=genoF=geno.1F=Popindex=NULL;
		uniPop=unique(ped[,'Population'])
		for(j in 1:length(uniPop))
		{
			tmp=geno[ped[,'Population']==uniPop[j],];
			Nfamtmp=floor(nrow(tmp)/(Nind));
			Ninctmp=Nfamtmp*Nind;
			tmp=tmp[sample(1:nrow(tmp),Ninctmp),];
			genoMtmp=tmp[1:Nfamtmp,]
			genoFtmp=tmp[(Nfamtmp+1):(2*Nfamtmp),]
			geno.1Ftmp=tmp[(2*Nfamtmp+1):nrow(tmp),];
			genoM=rbind(genoM,genoMtmp)
			genoF=rbind(genoF,genoFtmp)
			geno.1F=rbind(geno.1F,geno.1Ftmp);
			Popindex=rbind(Popindex,matrix(rep(uniPop[j],Ninctmp,ncol=1)))
		}

		geno.GF=genoF
		geno.GM=genoM
		Nfam=nrow(geno.GF);

	# The first generation #
		geno.1M=NULL;
		for(i in 1:Noff1)geno.1M=rbind(geno.1M,GeneDropCall(geno.GF,geno.GM));
	# The second generation #
		geno.2O=NULL;
		if(nrow(geno.1F)!=nrow(geno.1M)) {cat("Error","M!=F","\n");break;}
		for(i in 1:Noff2) geno.2O=rbind(geno.2O,GeneDropCall(geno.1F,geno.1M));
	# Building indexes #
		GENOSAVE=rbind(geno.GF,geno.GM,geno.1F,geno.1M,geno.2O);
		Father.id0=1:nrow(geno.GF);
		Mother.id0=(nrow(geno.GF)+1):nrow(rbind(geno.GF,geno.GM))
		Father.id1=(nrow(geno.GF)+nrow(geno.GM)+1):nrow(rbind(geno.GF,geno.GM,geno.1F))
		Mother.id1=(nrow(geno.GF)+nrow(geno.GM)+nrow(geno.1F)+1):nrow(rbind(geno.GF,geno.GM,geno.1F,geno.1M))

		IDSAVE=1:nrow(GENOSAVE);
		FATHERIDSAVE=c(rep(0,nrow(geno.GF)),rep(0,nrow(geno.GM)),rep(0,nrow(geno.1F)),rep(Father.id0,Noff1),rep(Father.id1,Noff2))
		MOTHERIDSAVE=c(rep(0,nrow(geno.GF)),rep(0,nrow(geno.GM)),rep(0,nrow(geno.1F)),rep(Mother.id0,Noff1),rep(Mother.id1,Noff2))
		FAMILYIDSAVE=c(1:Nfam,1:Nfam,rep(1:Nfam,Noff1),rep(1:Nfam,Noff1),rep(rep(1:Nfam,Noff1),Noff2));

		SEXSAVE=c(rep(1,nrow(geno.GF)),rep(2,nrow(geno.GM)),rep(1,nrow(geno.1F)),rep(2,nrow(geno.1M)),sample(c(1,2),nrow(geno.GF)*Noff1*Noff2,replace=T))
		PEDSAVE=cbind(IDSAVE,FATHERIDSAVE,MOTHERIDSAVE,SEXSAVE,FAMILYIDSAVE);
		colnames(PEDSAVE)=c("ID","DADID","MOMID","SEX","FAMID");
	}
	result=list();
	result$GENOSAVE=GENOSAVE;
	result$PEDSAVE=PEDSAVE;
	result;
}
	

