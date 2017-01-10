#library(SKAT)
library(CompQuadForm)
library(kinship2)

MAF.family<-function(Z,family)
{
	fam=unique(family[,'FAMID']);
	MAF.fam=matrix(-9,nrow=length(fam),ncol=ncol(Z));
	for(i in 1:length(fam))
	{
		t1=family[,'FAMID']==fam[i];
		geno.tmp=Z[t1,]
		if(sum(t1)==1) MAF.fam[i,]=geno.tmp/2
		if(sum(t1)>1) MAF.fam[i,]=apply(geno.tmp,2,mean,na.rm=T)/2;
	}
	colnames(MAF.fam)=colnames(Z);
	rownames(MAF.fam)=fam
	MAF.fam
}
null.FGRF.heter.pred<-function(Y,family,Z,genes,X=NULL,weight=NULL,out_type="C",MAF=NULL) #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix
{
	step.max=100;
	##Preliminary
	Y<-as.matrix(Y);n<-nrow(Y);index.cluster<-unique(family[,1])
	X1<-cbind(rep(1,n),X);p<-ncol(X1)-1;V<-diag(1,n);J<-matrix(0,n,n)
	par.geno<-Z[which(family[,3]<=0 & family[,4]<=0 ),];

	if(is.null(MAF)) MAF<-apply(par.geno,2,mean)/2 # calculate population based allele.freq;
	names(MAF)=colnames(Z);
	if(length(weight)==0){weight=rep(1,length(MAF))}#weight<-1/sqrt(2*MAF*(1-MAF))}#or dbeta(MAF,1,25) as in SKAT
	weight<-weight[MAF!=0];Z<-Z[,MAF!=0];MAF<-MAF[MAF!=0] # consider variants with variation

	# Family Specific allele Freq #
	MAF.fam=MAF.family(Z,family);



	#Get centered and weighted genotype
	#G<-Get.newG(Z,weight)
	S.genetic=list();
	S.genetic.fam=list();
	for(i in 1:length(genes))
	{
		S.genetic.fam[[i]]=matrix(0,n,n)
		S.genetic[[i]]=matrix(-9,n,n);
	}



	S.kinship<-matrix(0,n,n);S.neighbor<-matrix(0,n,n);S.uniform<-matrix(0,n,n);
	#within-family similarity
	Z.adj=matrix(-9,nrow(Z),ncol(Z));
	for (i in index.cluster)
	{
		index<-which(family[,1]==i);ni<-length(index);
		J[index,index]<-matrix(1/ni,ni,ni);
		S.kinship[index,index]<-as.matrix(makekinship(family[index,1], family[index,2], family[index,3], family[index,4]))
		S.neighbor[index,index]<-Distance.fmatrix(family[index,])
		S.uniform[index,index]<-matrix(1,ni,ni);
		
		for(j in 1:length(genes))
		{
			t1=col.index=colnames(Z) %in% genes[[j]]
			Z.tmp=Z[index,t1]
			sn=rank(colnames(Z.tmp))
			tt=colnames(MAF.fam) %in% colnames(Z.tmp);
			MAF.fam.tmp=MAF.fam[,tt];
			or=order(colnames(MAF.fam.tmp));
			MAF.fam.tmp=MAF.fam.tmp[,or][,sn];
			if(sum(colnames(MAF.fam.tmp)!=colnames(Z.tmp))!=0) print("ERROR");
			if(sum(colnames(MAF.fam.tmp)!=colnames(Z.tmp))==0) maf.fam.tmp=MAF.fam.tmp[rownames(MAF.fam.tmp)==i,];
			MAF.tmp=MAF[names(MAF) %in% colnames(Z.tmp)]
			or=order(names(MAF.tmp))
			maf.tmp=MAF.tmp[or][sn];
			
			Gi=t((t(Z.tmp)-maf.tmp*2)/sqrt(2*maf.tmp*(1-maf.tmp)));
			ss=Gi %*% t(Gi)
			z.adj=t((t(Z.tmp)-maf.fam.tmp*2)/sqrt(2*maf.tmp*(1-maf.tmp)));			
			Z.adj[index,col.index]=z.adj;
			S.genetic.fam[[j]][index,index]=ss;
			diag(S.genetic.fam[[j]][index,index])=0;
		}
		diag(S.kinship[index,index])<-0;diag(S.neighbor[index,index])<-0;diag(S.uniform[index,index])<-0
	}

	for(i in 1:length(genes))
	{
		t1=col.index=colnames(Z) %in% genes[[i]]
		Z.adj=Z.adj[,t1]
		S.genetic[[i]]=Z.adj %*% t(Z.adj);
		for (j in index.cluster)
		{
			index<-which(family[,1]==j)
			S.genetic[[i]][index,index]=0;
		#	diag(S.genetic[[i]])=0;
		}
	}

#	if((sum(unlist(S.genetic.fam)==(-9))+sum(unlist(S.genetic)==(-9))) !=0) print("Error")

	#Iteration to get beta and eta
	step.error<-Inf;beta<-rep(0,p+1);eta<-c(0,0,0);iter<-0 #iteration time
	eta<-c(0,0,rep(0,length(S.genetic)*2))

#	mu=rep(0,nrow(Y))
#	for (i in index.cluster)
#	{
#		index<-which(family[,1]==i);
#		mu[index]=Y[index]-mean(Y[index]);
#	}

#	print("RUN...");
	step.nn=0;
	step.beta.nn=0;
	while (step.error>0.00001 & step.nn<step.max)
	{
		step.nn=step.nn+1;
		step.beta.nn=0;
		step.beta<-Inf;
		mu<-rep(mean(Y),n);
		Y.res<-Y-mu;#initiate residual
		while (step.beta>0.00001 & step.beta.nn<step.max)
		{
			step.beta.nn=step.beta.nn+1;
			if (out_type=="C")
			{delta<-rep(1,n);A<-t(X1)%*%V%*%X1;B<-t(X1)%*%V%*%Y.res}
			if (out_type=="D")
			{
 				delta<-as.vector(mu*(1-mu));
				A<-t(delta^{1/2}*X1)%*%V%*%(delta^{1/2}*X1);
 				B<-t(delta^{1/2}*X1)%*%V%*%(delta^{-1/2}*Y.res)
			}
			#calculating residual
			beta.new<-beta+solve(A)%*%B;
			if (out_type=="C"){mu<-X1%*%beta.new;Y.res<-Y-mu}
			if (out_type=="D"){mu<-exp(X1%*%beta.new)/(1+exp(X1%*%beta.new));Y.res<-Y-mu}
			#calculate step beta error
			step.beta<-sum((beta.new-beta)^2);beta<-beta.new
		}
		#Generate RF psuedo variables
		#    x.field<-cbind(S.kinship%*%Y.res,S.neighbor%*%Y.res,S.uniform%*%Y.res)
		x.field<-cbind(S.kinship%*%Y.res,S.uniform%*%Y.res)
		x.field<-cbind(S.kinship%*%Y.res,S.uniform%*%Y.res)
		for(i in 1:length(genes))
		x.field<-cbind(x.field,S.genetic[[i]]%*%Y.res,S.genetic.fam[[i]]%*%Y.res)
		eta.new<-lm(Y.res~0+x.field)$coefficients
		#calculate step error
		step.error<-sum(c(beta.new-beta,eta.new-eta)^2);iter<-iter+1
		beta<-beta.new;eta<-eta.new
		#    V<-diag(1,n)-eta[1]*S.kinship-eta[2]*S.neighbor-eta[3]*S.uniform
		V<-diag(1,n)-eta[1]*S.kinship-eta[2]*S.uniform
		for(i in 1:length(genes))
		{
			V<-V-eta[(i-1)*2+2+1]*S.genetic[[i]]-eta[(i-1)*2+2+2]*S.genetic.fam[[i]]
		}
	}
	return(list(Converge=(step.nn<step.max),Y.res=Y.res,family=family,X1=X1,eta=eta,beta=beta,V=V,J=J,S.kinship=S.kinship,S.neighbor=S.neighbor,S.uniform=S.uniform,S.genetic=S.genetic,S.genetic.fam=S.genetic.fam,delta=delta,out_type=out_type))
}

