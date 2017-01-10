SampleStrategies<-function(PED,pfamily=2/3)
{
  FAMID=unique(PED[,"FAMID"]);
  sind=NULL;
  for(i in 1:length(FAMID))
  {
	tmp=PED[,"FAMID"]==FAMID[i]
	id=PED[tmp,"ID"];
	sind=c(sind,(sample(id,length(id)*pfamily)));
  }
  sind=sind[order(sind)]
  index=which(PED[,"ID"] %in% sind);
  index;
}
