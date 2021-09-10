

setwd('/Users/jamie/linux/diss/git/code/charts/VAR/')

require(readxl)
require(mFilter)

dt=read_excel('VARDATA.xlsx',sheet='Quarterly2')

ls=exp(dt$LS/100)

y=exp(dt$Y)


ls.hp=hpfilter(ls,type='lambda',freq=1600)

y.hp=hpfilter(log(y),type='lambda',freq=1600)

corr.y.ls=cor(ls.hp$cycle,y.hp$cycle)

print(corr.y.ls)

print(sd(ls))

print(sd(ls.hp$cycle))
