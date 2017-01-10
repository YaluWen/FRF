
set.seed(1)

peddir="./data/ped.infor"


progdir="./progfinal/";
chrsinput="./data/CHR1.RData";
source(paste(progdir,"FunctionSeqFam.r",sep=""))

nlocsave=3
nloctrue=1
effect=2
ndrop=1
#lprop=132
sd.err=1
pedscheme=2
Noff1=4
Noff2=0
pfamily=2/3
binary=FALSE
scheme="Equal"
out_type="C"

drop.effect=0.01

HeterAll=1/c(1,2,4,8,16,32,64,128)
#heter.prop=1
HeterMax=max(1/HeterAll);
wtinputs=c(0,1)

for(i in 1:length(wtinputs))
{
	Result=list();
	wtinput=wtinputs[i];
	cat(i,'\n')

	for(j in 1:50)
	{
		cat(j,'\n')
		for(k in 1:length(HeterAll))
		{
			heter.prop=HeterAll[k];lprop=HeterMax*heter.prop;
			if(k==1) 
			{
				GenoGenerate=FALSE;
				lprop=32;
				re=Simulation2(chrsinput=chrsinput,nlocsave=nlocsave,nloctrue=nloctrue,effect=effect,peddir=peddir,ndrop=ndrop,drop.effect=drop.effect,lprop=lprop,sd.err=sd.err,wtinput=wtinput,pedscheme=pedscheme,Noff1=Noff1,Noff2=Noff2,pfamily=pfamily,binary=binary,scheme=scheme,X=NULL,heter.prop=heter.prop)
				DataPast=list();
				DataPast$index=re$index;
				DataPast$model=re$model;
				if(j==1)	Result[[k]]=re$result;
				if(j!=1)	Result[[k]]=rbind(Result[[k]],re$result);
			}
			if(k!=1)
			{
				GenoGenerate=TRUE;
				if(k==2) lprop=16;
				re=Simulation2(chrsinput=chrsinput,nlocsave=nlocsave,nloctrue=nloctrue,effect=effect,peddir=peddir,ndrop=ndrop,drop.effect=drop.effect,lprop=lprop,sd.err=sd.err,wtinput=wtinput,pedscheme=pedscheme,Noff1=Noff1,Noff2=Noff2,pfamily=pfamily,binary=binary,scheme=scheme,X=NULL,heter.prop=heter.prop,GenoGenerate=GenoGenerate,DataPast=DataPast)
				if(j==1)	Result[[k]]=re$result;
				if(j!=1)	Result[[k]]=rbind(Result[[k]],re$result);
			}	
		}

	}

	for(k in 1:length(HeterAll))
	{
		write.table(Result[[k]],paste("He",1/HeterAll[k],"Drop",drop.effect,"Wt",wtinput,"Ped2",".txt",sep=""),quote=F,row.names=F)
	}	
}
