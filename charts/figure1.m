
loadd=1;

if loadd==1
clear all
close all

repdir='../../VAR/jvaa034_replication_codes_cantore_ferroni_leon_ledesma/replication_codes_Cantore_Ferroni_Leon_Ledesma/replication_codes_missing_link_appendix/Appendix A-B-D-F';

%load(strcat(repdir,'/Figure_B3/workspace.mat'))
load(strcat(repdir,'/Figures_B1_B8_B9_B10/workspace.mat'))

end

% US

% short term rate, log gdp, gdp def, CPI,PoC, W, LS

US_data=empirics_(1).dataest;
US_IRF=empirics_(1).bvar.ir_draws;
US_IRFproxy=empirics_(1).bvar.irproxy_draws_multiple;

meanLS=mean(exp(empirics_(1).dataest(:,7)/100));
normR=median(squeeze(US_IRFproxy(1,1,1,:))');


LSirf=100*(exp(squeeze(US_IRFproxy(7,:,1,:))'/normR/100)*meanLS-meanLS);

Wirf=squeeze(US_IRFproxy(6,:,1,:))'/normR;

YNirf=-squeeze(US_IRFproxy(7,:,1,:))'/normR+Wirf;

%% chart

figure(1)
clf

subplot(1,2,1)

lo=quantile(LSirf,0.5-0.68/2);
hi=quantile(LSirf,0.5+0.68/2);


   h = area([1:20],[lo',hi' - lo'],'FaceAlpha',0.20);
            
            set(h(2),'FaceColor',[0 0 .9])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
hold on
            plot([1:20],median(LSirf),'blue','LineWidth',1.8)


hline=refline(0,0);
hline.Color='black';
xlabel('Qtr')
ylabel('Percentage points')
title('Labour share')

xlim([1,20])

subplot(1,2,2)

lo=quantile(YNirf,0.5-0.68/2);
hi=quantile(YNirf,0.5+0.68/2);


   h = area([1:20],[lo',hi' - lo'],'FaceAlpha',0.20);
            
            set(h(2),'FaceColor',[0 0 .9])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
hold on
            plot([1:20],median(YNirf),'blue','LineWidth',1.8)


hline=refline(0,0);
hline.Color='black';
xlabel('Qtr')
ylabel('Percent')

title('Labour productivity')
xlim([1,20])

%% saving

figure(1)
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_US_LS','-dpdf','-r0')




