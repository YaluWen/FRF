#library(SKAT)
library(CompQuadForm)
library(kinship2)

#########################
#   Main function
#########################

#########################
#   Null model fitting
#########################

null.FGRF<-function(Y,family,X=NULL,out_type="C") #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix
{
	##Preliminary
	Y<-as.matrix(Y);n<-nrow(Y);index.cluster<-unique(family[,1])
	X1<-cbind(rep(1,n),X);p<-ncol(X1)-1;V<-diag(1,n);J<-matrix(0,n,n)
	S.kinship<-matrix(0,n,n);S.neighbor<-matrix(0,n,n);S.uniform<-matrix(0,n,n)
	
	#within-family similarity
	for (i in index.cluster)
	{
		index<-which(family[,1]==i);ni<-length(index)
		J[index,index]<-matrix(1/ni,ni,ni)
		S.kinship[index,index]<-as.matrix(makekinship(family[index,1], family[index,2], family[index,3], family[index,4]))
		S.neighbor[index,index]<-Distance.fmatrix(family[index,])
		S.uniform[index,index]<-matrix(1,ni,ni)
		diag(S.kinship[index,index])<-0;diag(S.neighbor[index,index])<-0;diag(S.uniform[index,index])<-0
	}

	#Iteration to get beta and eta
	step.error<-Inf;beta<-rep(0,p+1);eta<-c(0,0,0);iter<-0 #iteration time
	eta<-c(0,0);
	while (step.error>0.00001)
	{
		step.beta<-Inf;mu<-rep(mean(Y),n);Y.res<-Y-mu;#initiate residual
		while (step.beta>0.00001)
		{
			if (out_type=="C"){delta<-rep(1,n);A<-t(X1)%*%V%*%X1;B<-t(X1)%*%V%*%Y.res}
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
		eta.new<-lm(Y.res~0+x.field)$coefficients
		#calculate step error
		step.error<-sum(c(beta.new-beta,eta.new-eta)^2);iter<-iter+1
		beta<-beta.new;eta<-eta.new
		#    V<-diag(1,n)-eta[1]*S.kinship-eta[2]*S.neighbor-eta[3]*S.uniform
		V<-diag(1,n)-eta[1]*S.kinship-eta[2]*S.uniform
	}	
	return(list(Y.res=Y.res,family=family,X1=X1,eta=eta,beta=beta,V=V,J=J,S.kinship=S.kinship,S.neighbor=S.neighbor,S.uniform=S.uniform,delta=delta,out_type=out_type))
}

adjust.pred<-function(Y.tr,family.tr,Z.tr,Z.te,Z,index,genes,X.tr=NULL,X.te=NULL,out_type="C",weight=1)
{
	null=null.FGRF(Y=Y.tr,family=family.tr,X=X.tr,out_type=out_type)
	x.field<-cbind(null$S.kinship%*%null$Y.res,null$S.uniform%*%null$Y.res)
	y.null=x.field %*% matrix(null$eta,ncol=1)
	y.gene=null$Y.res-y.null;
	S.genetic=S.genetic.tr=S.genetic.te=list();
	x.field.genetic.tr=NULL;
	x.field.genetic.te=NULL;

	for(i in 1:length(genes))
	{
		S.genetic[[i]]=getIBS(Z,weight)
		S.genetic.tr[[i]]=S.genetic[[i]][index,index]
		S.genetic.te[[i]]=S.genetic[[i]][-index,index];
		x.field.genetic.tr=cbind(x.field.genetic.tr,S.genetic.tr[[i]] %*% y.gene)
		x.field.genetic.te=cbind(x.field.genetic.te,S.genetic.te[[i]] %*% y.gene)

	}
	pre.effect=lm(y.gene~x.field.genetic.tr-1)$coef
	pred=x.field.genetic.te %*% matrix(pre.effect,ncol=1)
	pred
}

#######################
##### My method   #####
#######################
null.FGRF.pred<-function(Y,family,Z,genes,X=NULL,weight=NULL,out_type="C") #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix
{
	##Preliminary
	Y<-as.matrix(Y);n<-nrow(Y);index.cluster<-unique(family[,1])
	X1<-cbind(rep(1,n),X);p<-ncol(X1)-1;V<-diag(1,n);J<-matrix(0,n,n)
	par.geno<-Z[which(family[,3]<=0 & family[,4]<=0 ),];

	MAF<-apply(par.geno,2,mean)/2 # calculate population based allele.freq
	if(length(weight)==0){weight=rep(1,length(MAF))}#weight<-1/sqrt(2*MAF*(1-MAF))}#or dbeta(MAF,1,25) as in SKAT
	weight<-weight[MAF!=0];Z<-Z[,MAF!=0];MAF<-MAF[MAF!=0] # consider variants with variation

	#Get centered and weighted genotype
	#G<-Get.newG(Z,weight)
	S.genetic=list()
	for(i in 1:length(genes))
	{
############################################### MODIFICATION ############
		col.index=colnames(Z) %in% genes[[i]]
		S.genetic[[i]]=getIBS(Z[,col.index],weight[col.index])
	}

	S.kinship<-matrix(0,n,n);S.neighbor<-matrix(0,n,n);S.uniform<-matrix(0,n,n);
	#within-family similarity
	for (i in index.cluster)
	{
		index<-which(family[,1]==i);ni<-length(index)
		J[index,index]<-matrix(1/ni,ni,ni)
		S.kinship[index,index]<-as.matrix(makekinship(family[index,1], family[index,2], family[index,3], family[index,4]))
		S.neighbor[index,index]<-Distance.fmatrix(family[index,])
		S.uniform[index,index]<-matrix(1,ni,ni)
		diag(S.kinship[index,index])<-0;diag(S.neighbor[index,index])<-0;diag(S.uniform[index,index])<-0
	}
	#Iteration to get beta and eta
	step.error<-Inf;beta<-rep(0,p+1);eta<-c(0,0,0);iter<-0 #iteration time
	eta<-c(0,0,rep(0,length(S.genetic)))
	while (step.error>0.00001)
	{
		step.beta<-Inf;mu<-rep(mean(Y),n);Y.res<-Y-mu;#initiate residual
		while (step.beta>0.00001)
		{
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
		x.field<-cbind(x.field,S.genetic[[i]]%*%Y.res)
		eta.new<-lm(Y.res~0+x.field)$coefficients
		#calculate step error
		step.error<-sum(c(beta.new-beta,eta.new-eta)^2);iter<-iter+1
		beta<-beta.new;eta<-eta.new
		#    V<-diag(1,n)-eta[1]*S.kinship-eta[2]*S.neighbor-eta[3]*S.uniform
		V<-diag(1,n)-eta[1]*S.kinship-eta[2]*S.uniform
		for(i in 1:length(genes))
		{
			V<-V-eta[i+2]*S.genetic[[i]]
		}
	}
	return(list(Y.res=Y.res,family=family,X1=X1,eta=eta,beta=beta,V=V,J=J,S.kinship=S.kinship,S.neighbor=S.neighbor,S.uniform=S.uniform,S.genetic=S.genetic,delta=delta,out_type=out_type))
}


############################################################
################  Other required function   ################
############################################################
Distance.fmatrix<-function(X) #family distance matrix, with diag = 0
{
  #checking parent-kid
  temp<-outer(X[,2],X[,3],Check.share)+outer(X[,2],X[,4],Check.share)+outer(X[,3],X[,2],Check.share)+outer(X[,4],X[,2],Check.share)
  return(temp)
}

Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){
    stop("Error to obtain matrix square!")
  }
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}

Get.newG<-function(Z,weights) #get centered and weighted genotype
{
  Make.center<-function(x){return(x-mean(x))}
  G<-apply(Z,2,Make.center)%*%diag(weights)
  return(G)
}
