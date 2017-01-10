
Simulation1<-function(chrsinput=chrsinput,nlocsave=nlocsave,nloctrue=nloctrue,effect=effect,ndrop=ndrop,drop.effect=drop.effect,
lprop=lprop,sd.err=sd.err,wtinput=wtinput,pedscheme=pedscheme,Noff1=Noff1,Noff2=Noff2,pfamily=pfamily,binary=binary,scheme=scheme,X=NULL,out_type="C")
{
	model=DiscodeFinal(chrsinput=chrsinput,nlocsave=nlocsave,nloctrue=nloctrue,effect=effect,ndrop=ndrop,drop.effect=drop.effect,lprop=lprop,sd.err=sd.err,wtinput=wtinput,pedscheme=pedscheme,Noff1=Noff1,Noff2=Noff2,pfamily=pfamily,binary=binary,scheme=scheme)
	ped=model$PEDSAVE;
	Z=model$GENOSAVE;
	Y=model$DIS;
	genes=Gene=model$GENE;

	family=cbind(ped[,'FAMID'],ped[,'ID'],ped[,'DADID'],ped[,'MOMID']);
	colnames(family)=cbind('FAMID','ID','DADID','MOMID');	

	index=SampleStrategies(PED=model$PEDSAVE,pfamily=pfamily)
	Y.tr=matrix(Y[index],ncol=1)
	Y.te=matrix(Y[-index],ncol=1)
	family.tr=family[index,]
	Z.tr=Z[index,]
	Z.te=Z[-index,]

	X=NULL
	X.te=X.tr=NULL;
	weight=1;
	AUCTR=AUCTE=AUCAD=NULL;
	#######################################################
	###################### MY METHOD  ######################
	########################################################

	result.null.pred=null.FGRF.pred(Y,family,Z,genes,X,weight=NULL,out_type=out_type) #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix		
	result.null.pred.tr=null.FGRF.pred(Y.tr,family.tr,Z.tr,genes,X,weight=NULL,out_type=out_type) #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix
	aa=matrix(cbind(rep(1,nrow(Y.tr)),X),ncol=1)
	Y.res=Y.tr-mean(Y.tr)
	aa=cbind(aa,result.null.pred.tr$S.kinship %*% Y.res,result.null.pred.tr$S.uniform%*% Y.res)
	for(i in 1:length(genes))
	aa=cbind(aa,result.null.pred.tr$S.genetic[[i]]%*% Y.res)

	bb=c(result.null.pred.tr$beta,result.null.pred.tr$eta)
	bb=matrix(bb,ncol=1)
	Ypred=aa %*% bb
	if(!binary)	auc=cor(Ypred,Y.tr)
	if(binary)
	{
		auc=wilcox.test(Ypred[Y.tr==0],Ypred[Y.tr==1])$statistic/sum(Y.tr)/sum(Y.tr==0)
		auc=max(auc,1-auc)
	}
	AUCTR=rbind(AUCTR,auc);
	
	aa=matrix(cbind(rep(1,nrow(Y.te)),X),ncol=1)
	Y.res=Y.tr-mean(Y.tr)
	aa=cbind(aa,result.null.pred$S.kinship[-index,index] %*% Y.res,result.null.pred$S.uniform[-index,index]%*% Y.res)
	for(i in 1:length(genes))
	aa=cbind(aa,result.null.pred$S.genetic[[i]][-index,index]%*% Y.res)
	bb=c(result.null.pred.tr$beta,result.null.pred.tr$eta)
	bb=matrix(bb,ncol=1)
	Ypred=aa %*% bb
	if(binary)
	{
		auc=wilcox.test(Ypred[Y.te==0],Ypred[Y.te==1])$statistic/sum(Y.te)/sum(Y.te==0)
		auc=max(auc,1-auc)
	}
	if(!binary) auc=cor(Ypred,Y.te)
	AUCTE=rbind(AUCTE,auc);

	#######################################################
	################# Adjusted Ones #######################
	#######################################################

	Ypred1=adjust.pred(Y.tr,family.tr,Z.tr,Z.te,Z,index,genes,X.tr=X.tr,X.te=X.te,out_type=out_type,weight=weight)
	if(binary)
	{
		auc=wilcox.test(Ypred1[Y.te==0],Ypred1[Y.te==1])$statistic/sum(Y.te)/sum(Y.te==0)
		auc=max(auc,1-auc)
	}
	if(!binary) auc=cor(Ypred1,Y.te)
	AUCAD=rbind(AUCAD,auc);

	result=c(AUCTR,AUCTE,AUCAD);
	names(result)=c("TR","TE","AD")
	result
}
	
	
