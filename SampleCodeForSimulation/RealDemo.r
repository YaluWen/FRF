progdir="./progfinal/";
source(paste(progdir,"FunctionSeqFam.r",sep=""))
## genes is a list with each element represent a chunk of genome (e.g. gene)
## family is the pedigree informaion, with family ID, ID, dad ID, and mom ID.
## Y is the outcome
## Z is a matrix containing all genotypes 
## index representing which data is used to train the model

# Y.tr=matrix(Y[index],ncol=1) , outcomes for training data set
# Y.te=matrix(Y[-index],ncol=1), outcomes for testing data 
# family.tr=family[index,], family informatino for training data
# Z.tr=Z[index,], genotypes for training data
# Z.te=Z[-index,], genotypes for testing data
# out_type="C", continuous outcomes, and out_type="D" is for binary outcomes #

# X=NULL, other covariates #,
# X.te=X.tr=NULL, other covariates #
# weight=1, weight used to reflect the contribution of rare variants, by default it is equal weight #


## Using all data to fit the model ##
result.null.pred=null.FGRF.pred(Y,family,Z,genes,X,weight=NULL,out_type=out_type) #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix		

## Use only training data to fit the model ##
result.null.pred.tr=null.FGRF.pred(Y.tr,family.tr,Z.tr,genes,X,weight=NULL,out_type=out_type) #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix
aa=matrix(cbind(rep(1,nrow(Y.tr)),X),ncol=1)
Y.res=Y.tr-mean(Y.tr)
aa=cbind(aa,result.null.pred.tr$S.kinship %*% Y.res,result.null.pred.tr$S.uniform%*% Y.res)
for(i in 1:length(genes))  aa=cbind(aa,result.null.pred.tr$S.genetic[[i]]%*% Y.res)
bb=c(result.null.pred.tr$beta,result.null.pred.tr$eta)
bb=matrix(bb,ncol=1)
Ypred=aa %*% bb  ## predicted value for training data 
  
aa=matrix(cbind(rep(1,nrow(Y.te)),X),ncol=1)
Y.res=Y.tr-mean(Y.tr)
aa=cbind(aa,result.null.pred$S.kinship[-index,index] %*% Y.res,result.null.pred$S.uniform[-index,index]%*% Y.res)
for(i in 1:length(genes)) aa=cbind(aa,result.null.pred$S.genetic[[i]][-index,index]%*% Y.res)
bb=c(result.null.pred.tr$beta,result.null.pred.tr$eta)
bb=matrix(bb,ncol=1)
Ypred=aa %*% bb ## predicted values for testing data using the model calculated from training data
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


