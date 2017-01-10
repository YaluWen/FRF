
library(Rcpp)
library(inline)
# genoF: The genotype of the father with its row representing individual, and column represents the SNP #
# genoM: The genotype of the mother with its row representing individual, and column represents the SNP #
src=
  '
NumericVector Seed(seed);
NumericMatrix GenoF(genoF);
NumericMatrix GenoM(genoM);
NumericMatrix GenoO(genoO);

if(Seed[0]==1)
{
	int N=Seed[0];
	srand(N);
}

int nrowF=GenoF.nrow();
int ncolF=GenoF.ncol();
int nrowM=GenoM.nrow();
int ncolM=GenoM.ncol();

int run=1;
if(nrowF!=nrowM){std::cout<<"Number of Father Does not Equal to Number oF Mother"<<std::endl;run=0;}
if(ncolF!=ncolM){std::cout<<"Number of Father\'s geno Does not Equal to Number oF Mother\'s geno"<<std::endl;run=0;}

if(run==1)
{
	NumericVector fathertmp(ncolF,0.0);
	NumericVector mothertmp(ncolM,0.0);
	int tmp=0;
	for(int i=0; i<nrowF;i++)
	{
		fathertmp=GenoF(i,_);
		mothertmp=GenoM(i,_);
		for(int j=0;j<ncolF;j++)
		{
			GenoO(i,j)=-9;
			if(fathertmp[j]==0 && mothertmp[j]==0) GenoO(i,j)=0;
			if(fathertmp[j]==0 && mothertmp[j]==1) GenoO(i,j)=rand()%2;
			if(fathertmp[j]==0 && mothertmp[j]==2) GenoO(i,j)=1;
			if(fathertmp[j]==1 && mothertmp[j]==0) GenoO(i,j)=rand()%2;
			if(fathertmp[j]==1 && mothertmp[j]==1) 
			{
				tmp=rand()%4;
				GenoO(i,j)=1;
				if(tmp==0) GenoO(i,j)=0;
				if(tmp==3) GenoO(i,j)=2;
			}
			if(fathertmp[j]==1 && mothertmp[j]==2) GenoO(i,j)=rand()%2+1;
			if(fathertmp[j]==2 && mothertmp[j]==0) GenoO(i,j)=1;
			if(fathertmp[j]==2 && mothertmp[j]==1) GenoO(i,j)=rand()%2+1;
			if(fathertmp[j]==2 && mothertmp[j]==2) GenoO(i,j)=2;
		}
	}
}

return Rcpp::wrap(GenoO);
'	

GeneDrop=cxxfunction(signature(genoF='matrix',genoM='matrix',genoO='matrix',seed='numeric'),src,plugin="Rcpp")

genoF=matrix(c(0,1,2,1,2,0),nrow=2);
genoM=matrix(c(0,12,2,1,1,2),nrow=2);

GeneDropCall<-function(genoF,genoM,seed=sample(1:10000,1))
{
	genoO=matrix(,nrow=nrow(genoM),ncol=ncol(genoM))
	genoO=GeneDrop(genoF,genoM,genoO,seed)
	genoO
}
