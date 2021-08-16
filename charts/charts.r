

setwd('/media/jamie/OS/Users/jamie/linux/diss/git/MResDiss/')

require(ggplot2)
require(tidyr)
require(lubridate)
require(data.table)
require(viridis)
require(gridExtra)

dt=data.table(read.csv('decompdata.csv',header = FALSE))

colnames(dt)=c('Qtr','HANK','HANKYN')

dt$Var=c(rep('Total',40),rep('rb',40),rep('ra',40),rep('Lab',40),rep('Tax',40))

fga=ggplot(data=dt[Var!='Total'])+geom_bar(aes(x=Qtr,y=HANK,fill=factor(Var)),stat='identity')+
  geom_line(data=dt[Var=='Total'],aes(x=Qtr,y=HANK,color='Total'),size=1.2,linetype=2)+
  labs(y='PP',fill='',color='',title='HANK')+
  theme_classic()+scale_colour_manual(values=c('black'))+scale_fill_viridis_d()+
  guides(fill="none",color="none")+scale_y_continuous(breaks=seq(-2,2,0.2),limits = c(-1.2,0.4))+
  geom_hline(yintercept = 0)

fgb=ggplot(data=dt[Var!='Total'])+geom_bar(aes(x=Qtr,y=HANKYN,fill=factor(Var)),stat='identity')+
  geom_line(data=dt[Var=='Total'],aes(x=Qtr,y=HANKYN,color='Total'),size=1.2,linetype=2)+
  labs(y='',fill='',color='',title='HANK-YN')+
  theme_classic()+scale_colour_manual(values=c('black'))+scale_fill_viridis_d()+
  theme(legend.position = c(0.75,0.30),legend.key = element_rect(colour = "transparent", fill = "white"))+
  scale_y_continuous(breaks=seq(-2,2,0.2),limits = c(-1.2,0.4))+
  geom_hline(yintercept = 0)

fgdecomp=grid.arrange(fga,fgb,nrow=1)

ggsave('decomp.pdf',fgdecomp,units = 'cm',width=16,height = 12)
  