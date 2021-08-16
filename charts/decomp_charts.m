clc
clear all
close all

addpath(genpath('../HANKYN/functions'))

%% Decompose NK consumption

%partial_decomp_NK(SS,T,Tf,Irb,Ira,Ilab,Itau)

[IRF_C_NKtot,IRF_C]=partial_decomp_NK('../../../steadystates/NKfund_60_15.mat',101,40,1,1,1,1);
[IRF_C_NKrb]=partial_decomp_NK('../../../steadystates/NKfund_60_15.mat',101,40,1,0,0,0);
[IRF_C_NKra]=partial_decomp_NK('../../../steadystates/NKfund_60_15.mat',101,40,0,1,0,0);
[IRF_C_NKlab]=partial_decomp_NK('../../../steadystates/NKfund_60_15.mat',101,40,0,0,1,0);
[IRF_C_NKtax]=partial_decomp_NK('../../../steadystates/NKfund_60_15.mat',101,40,0,0,0,1);

scal1=IRF_C(1)./IRF_C_NKtot(1);
%scal=ones(1,24);

figure(1)
clf

subplot(1,2,1)

hold on 
plot(IRF_C_NKtot*scal1,'LineWidth',1.8)
plot(IRF_C_NKrb*scal1,'LineWidth',1.8)
plot(IRF_C_NKra*scal1,'LineWidth',1.8)
plot(IRF_C_NKlab*scal1,'LineWidth',1.8)
plot(IRF_C_NKtax*scal1,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

legend('Tot','rb','ra','Labour','Tax','Location','southeast')

title('HANK')

%% Decompse YN consumption
 

[IRF_C_YNtot,IRF_C]=partial_decomp_YN('../../../steadystates/YNfund_60_15.mat',101,40,1,1,1,1);
[IRF_C_YNrb]=partial_decomp_YN('../../../steadystates/YNfund_60_15.mat',101,40,1,0,0,0);
[IRF_C_YNra]=partial_decomp_YN('../../../steadystates/YNfund_60_15.mat',101,40,0,1,0,0);
[IRF_C_YNlab]=partial_decomp_YN('../../../steadystates/YNfund_60_15.mat',101,40,0,0,1,0);
[IRF_C_YNtax]=partial_decomp_YN('../../../steadystates/YNfund_60_15.mat',101,40,0,0,0,1);

scal2=IRF_C(1)/IRF_C_YNtot(1);
%scal=ones(1,24);

figure(1)


subplot(1,2,2)
hold on 
plot(IRF_C_YNtot*scal2,'LineWidth',1.8)
plot(IRF_C_YNrb*scal2,'LineWidth',1.8)
plot(IRF_C_YNra*scal2,'LineWidth',1.8)
plot(IRF_C_YNlab*scal2,'LineWidth',1.8)
plot(IRF_C_YNtax*scal2,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

legend('Tot','rb','ra','Labour','Tax','Location','southeast')


title('HANK-YN')

%%
dtout=[ repmat( [1:40]',5,1), scal1*[IRF_C_NKtot'; IRF_C_NKrb'; IRF_C_NKra'; IRF_C_NKlab'; IRF_C_NKtax'], scal2*[IRF_C_YNtot'; IRF_C_YNrb'; IRF_C_YNra'; IRF_C_YNlab'; IRF_C_YNtax']];


csvwrite('decompdata.csv',dtout)


