clear all

close all

cd('../../../RepAgent/estimate/')

load('DynamicsNK_0827_IRF.mat')

NKcol =  [0.2 0.2 1.0];
YNcol =  [1 0 0.83];

figure(1)
clf
scalNK=0.01/oo_.posterior_mean.shocks_std.eint;
hold on
h1=plot(100*oo_.PosteriorIRF.dsge.Median.sy_eint*scalNK,'Color',NKcol,'LineWidth',1.8);


loNK=100*oo_.PosteriorIRF.dsge.HPDinf.sy_eint*scalNK;
hiNK=100*oo_.PosteriorIRF.dsge.HPDsup.sy_eint*scalNK;

hline=refline(0,0);
hline.Color='black';
ylabel('Percentage points')

%%
figure(2)
clf
subplot(1,2,1)
hold on
scal=1;
%scal=oo_.posterior_mean.shocks_std.eint;
plot(100*oo_.PosteriorIRF.dsge.Median.sy_eint/scal,'LineStyle','-','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ea;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_ezy/scal,'LineStyle','--','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezi;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_ezi/scal,'LineStyle','--','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.egov;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_egov/scal,'LineStyle','-','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezrp;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_ezrp/scal,'LineStyle','-','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezmk;
plot(100*oo_.PosteriorIRF.dsge.Median.sy_ezmk/scal,'LineStyle','--','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezwmk;
plot(100*oo_.PosteriorIRF.dsge.Median.sy_ezwmk/scal,'LineStyle','--','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
ylabel('Percentage points')
title('NK')
txt = 'Pro-cyclical';
text(20,-0.2,txt)


txt = 'Counter-cyclical';
text(20,0.2,txt)



%%
load('DynamicsNKYN_0827_IRF.mat')

figure(1)

scalYN=0.01/oo_.posterior_mean.shocks_std.eint;
hold on
h2=plot(100*oo_.PosteriorIRF.dsge.Median.sy_eint*scalNK,'Color',YNcol,'LineWidth',1.8);

      legend('NK','NK-YN','Location','southeast')

   hNK = area([1:40],[loNK,hiNK - loNK],'FaceAlpha',0.20);
            
            set(hNK(2),'FaceColor',NKcol)
            set(hNK(1),'FaceColor',[1 1 1])
            set(hNK,'linestyle','none')
            
lo=100*oo_.PosteriorIRF.dsge.HPDinf.sy_eint*scalNK;
hi=100*oo_.PosteriorIRF.dsge.HPDsup.sy_eint*scalNK;

   h = area([1:40],[lo,hi - lo],'FaceAlpha',0.20);
            
            set(h(2),'FaceColor',YNcol)
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')

legend([h1 h2],{'NK','NK-YN'})
figure(2)

subplot(1,2,2)
hold on
scal=1;
%scal=oo_.posterior_mean.shocks_std.eint;
plot(100*oo_.PosteriorIRF.dsge.Median.sy_eint/scal,'LineStyle','-','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezy;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_ezy/scal,'LineStyle','--','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezi;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_ezi/scal,'LineStyle','--','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.egov;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_egov/scal,'LineStyle','-','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezrp;
plot(-100*oo_.PosteriorIRF.dsge.Median.sy_ezrp/scal,'LineStyle','-','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezmk;
plot(100*oo_.PosteriorIRF.dsge.Median.sy_ezmk/scal,'LineStyle','--','LineWidth',1.8)
%scal=oo_.posterior_mean.shocks_std.ezwmk;
plot(100*oo_.PosteriorIRF.dsge.Median.sy_ezwmk/scal,'LineStyle','--','LineWidth',1.8)


hline=refline(0,0);
hline.Color='black';

legend('MP','TFP','IST','Gov','RP','Markup','Wage')
title('NK-YN')


figure(3)
clf
hold on
   h = area(oo_.posterior_density.parameters.shareE(:,1),oo_.posterior_density.parameters.shareE(:,2),'FaceAlpha',0.20);
   h2 = area(oo_.prior_density.parameters.shareE(:,1),oo_.prior_density.parameters.shareE(:,2),'FaceAlpha',0.20);
             set(h2,'FaceColor',[0.2 0.7 0])
             set(h,'linestyle','--','facealpha',0.75)
             %set(h2,'linestyle','none')
legend('Posterior','Prior')
   hold on
   xlim([0 1])
%xline(oo_.posterior_mean.parameters.shareE);

%% Tables


load('DynamicsNK_0827_VAR.mat')


vdyNK=   oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean.dy;

pNK=[struct2array(oo_.posterior_hpdinf.parameters)';struct2array(oo_.posterior_hpdinf.shocks_std)';
    struct2array(oo_.posterior_mean.parameters)';struct2array(oo_.posterior_mean.shocks_std)';
    struct2array(oo_.posterior_mode.parameters)';struct2array(oo_.posterior_mode.shocks_std)';
    struct2array(oo_.posterior_hpdsup.parameters)';struct2array(oo_.posterior_hpdsup.shocks_std)';];

ss=size(struct2array(oo_.posterior_hpdinf.parameters),2)+7;

pNK=reshape(pNK,ss,4);

fnNK=[fieldnames(oo_.posterior_mean.parameters);fieldnames(oo_.posterior_mean.shocks_std)];
pNK=table(fnNK,pNK);

%csvwrite('paramestNK.csv',pNK)

load('DynamicsNKYN_0827_VAR.mat')

vdyNKYN=   oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean.dy;

pNKYN=[struct2array(oo_.posterior_hpdinf.parameters)';struct2array(oo_.posterior_hpdinf.shocks_std)';
    struct2array(oo_.posterior_mean.parameters)';struct2array(oo_.posterior_mean.shocks_std)';
    struct2array(oo_.posterior_mode.parameters)';struct2array(oo_.posterior_mode.shocks_std)';
    struct2array(oo_.posterior_hpdsup.parameters)';struct2array(oo_.posterior_hpdsup.shocks_std)';];

ss=size(struct2array(oo_.posterior_hpdinf.parameters),2)+7;

pNKYN=reshape(pNKYN,ss,4);

fnNKYN=[fieldnames(oo_.posterior_mean.parameters);fieldnames(oo_.posterior_mean.shocks_std)];
pNKYN=table(fnNKYN,pNKYN);

%csvwrite('paramestNKYN.csv',pNKYN)
%% save figures

cd('../../git/code/charts/')

figure(1)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
%h.Position(3) =h.Position(3)*1.2;
%h.Position(4) =h.Position(4)*0.9;
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_lsest','-dpdf','-r0')

figure(2)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
%h.Position(3) =h.Position(3)*1.2;
%h.Position(4) =h.Position(4)*0.9;
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'irf_allshocks','-dpdf','-r0')

figure(3)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
%h.Position(3) =h.Position(3)*1.2;
%h.Position(4) =h.Position(4)*0.9;
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'shareE','-dpdf','-r0')
