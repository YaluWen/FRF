set.seed(10)


progdir="./progfinal/";
chrsinput="./data/CHR1.RData";
source(paste(progdir,"FunctionSeqFam.r",sep=""))


nlocsave=2
nloctrue=1
effect=0.2
ndrop=1
lprop=3
sd.err=1
pedscheme=2
Noff1=2
Noff2=0
pfamily=2/3
binary=FALSE
scheme="Equal"
out_type="C"

drop.effect=0.2


wtinputs=c(0,1,2,4,16,30)
Result=NULL;
for(i in 1:length(wtinputs))
{
  wtinput=wtinputs[i];
  cat(i,'\n')
  for(j in 1:1000)
    Result=rbind(Result,Simulation1(chrsinput=chrsinput,nlocsave=nlocsave,nloctrue=nloctrue,effect=effect,ndrop=ndrop,drop.effect=drop.effect,lprop=lprop,sd.err=sd.err,wtinput=wtinput,pedscheme=pedscheme,Noff1=Noff1,Noff2=Noff2,pfamily=pfamily,binary=binary,scheme=scheme,X=NULL))
  write.table(Result,paste("Drop",drop.effect,"Wt",wtinput,"Ped2",".txt",sep=""),quote=F,row.names=F)
  Result=NULL;
}
