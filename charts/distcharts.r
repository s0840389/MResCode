
require(ggplot2)
require(data.table)
require(tidyr)
require(viridis)

setwd('/media/jamie/OS/Users/jamie/linux/diss/git/MResDiss/')

readd=1

if(readd==1){
dt=read.csv('ynfund_dist.csv',header = FALSE)
dt=data.table(dt)
}

colnames(dt)=c('wgt','wealth','k','m','h','Wpct',paste('cirf',c(1:39),sep='_'),paste('clevelirf',c(1:39),sep='_'),paste('yirf',c(1:39),sep='_'))


fg=ggplot(data=dt)+geom_point(aes(x=k,y=cirf_1,color=m))+scale_color_viridis()


dtl=melt(dt,id.var=c('wgt','wealth','k','m','h','Wpct'),measure.vars = patterns('cirf_','clevelirf_','yirf',cols=names(dt)))


dtl$t=as.numeric(dtl$variable)

fg=ggplot(data=dtl)+geom_line(aes(x=t,y=value1,color=Wpct))+scale_color_viridis()
