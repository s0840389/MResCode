

clear all
close all

NKcol =  [0.2 0.2 1.0];
YNcol =  [1 0 0.83];
YN2Wcol= [0.0 0.8 0.0];

%% IRF-HANK

% NK model
load('../../steadystates/NKfund_60_15.mat')

IRF_NK=IRF_state_sparse;

figure(3)
clf

subplot(4,2,1)

plot(IRF_Y,'Color',NKcol,'LineWidth',1.8)

subplot(4,2,2)

plot(IRF_PI/100,'Color',NKcol,'LineWidth',1.8)

subplot(4,2,3)

plot(IRF_C,'Color',NKcol,'LineWidth',1.8)

subplot(4,2,4)

plot(IRF_I,'Color',NKcol,'LineWidth',1.8)

subplot(4,2,5)

plot(IRF_N,'Color',NKcol,'LineWidth',1.8)

subplot(4,2,6)

plot(IRF_W,'Color',NKcol,'LineWidth',1.8)

subplot(4,2,7)

plot(100*(IRF_LS-par.Labshr),'Color',NKcol,'LineWidth',1.8)

subplot(4,2,8)

plot(IRF_Y-IRF_N,'Color',NKcol,'LineWidth',1.8)

%%

figure(4) % otehr variables

clf

subplot(1,3,1)
plot(0.25*targets.BY*(IRF_B-IRF_Y),'Color',NKcol,'LineWidth',1.8)

subplot(1,3,2)
plot(IRF_A,'Color',NKcol,'LineWidth',1.8)

subplot(1,3,3)
hold on
plot(IRF_ra,'Color',NKcol,'LineWidth',1.8)


%%
figure(6)
clf
newcolors = [0.99 0.0 0.2;
    0.90 0.00 0.35;
    0.80 0.00 0.50;
    0.70 0.05 0.65;
    0.60 0.00 0.8;
    0.5 0.0 0.95];
           
colororder(newcolors)
            
subplot(1,2,1)

plot(IRF_hhc(1:6,:)'*100,'LineWidth',1.8)
hold on
hline=refline(0,0);
hline.Color='black';

title('HANK')


%%
figure(7)
clf

plot(IRF_giniC,'Color',NKcol,'LineWidth',1.8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% YN model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../../steadystates/YNfund_60_15.mat')

IRF_YN=IRF_state_sparse;

IRF_W_YN=IRF_W;
IRF_LS_YN=IRF_LS;
IRF_Y_YN=IRF_Y;
IRF_N_YN=IRF_N;
mpar_YN=mpar;
figure(3)

%output
subplot(4,2,1)
hold on

plot(IRF_Y,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Output gap')

legend('NK','NK-YN','Location','southeast')

ylabel('PP')

%inflation
subplot(4,2,2)

hold on

plot(IRF_PI/100,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Inflation')

%consumption
subplot(4,2,3)

hold on
plot(IRF_C,'Color',YNcol,'LineWidth',1.8)


hline=refline(0,0);
hline.Color='black';
title('Consumption')

ylabel('PP')

% investment
subplot(4,2,4)

hold on
plot(IRF_I,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Investment')

%Labour
subplot(4,2,5)

hold on

plot(IRF_N,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Hours')
ylabel('PP')

%Real wages
subplot(4,2,6)

hold on

plot(IRF_W,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Real wages')

%Labour share
subplot(4,2,7)

hold on
plot((IRF_LS-par.Labshr)*100,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Labour share')
ylabel('PP')

%Productivity
subplot(4,2,8)

hold on

plot(IRF_Y-IRF_N,'Color',YNcol,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Labour productivity')

%%

figure(4) % other variables

subplot(1,3,1)

hold on
plot(0.25*targets.BY*(IRF_B-IRF_Y),'Color',YNcol,'LineWidth',1.8)

subplot(1,3,2)

hold on
plot(IRF_A,'Color',YNcol,'LineWidth',1.8)

subplot(1,3,3)

hold on
plot(IRF_ra,'Color',YNcol,'LineWidth',1.8)


%%

figure(6)

subplot(1,2,2)
plot(IRF_hhc(1:6,:)'*100,'LineWidth',1.8)
hold on

hline=refline(0,0);
hline.Color='black';

legend(string(round(targets.wpcts*100,0)),'Location','southeast')
title('HANK-YN')

%% 

figure(7)

hold on
plot(IRF_giniC,'Color',YNcol,'LineWidth',1.8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IRF-HANK2w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('../../steadystates/YNfund2W_60_15.mat')

IRF_YN2W=IRF_state_sparse;

figure(5)
clf

newcolors = [YN2Wcol;
             YNcol];

            colororder(newcolors)
            
            
subplot(4,2,1)


hold on

title('Output gap')

plot(IRF_Y,'LineWidth',1.8)

plot(IRF_YN(end-mpar_YN.oc+3,1:end-1)*100,'','LineWidth',1.8)

legend('$w_e \ne w_y$','$w_e = w_y$','Location','southeast','AutoUpdate','off','Interpreter','latex')

hline=refline(0,0);
hline.Color='black';

ylabel('PP')

subplot(4,2,2)


plot(IRF_PI/100,'LineWidth',1.8)

hold on

title('Inflation')

plot(IRF_YN(end-mpar_YN.oc+2,1:end-1)*100,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

subplot(4,2,3)

hold on
plot(IRF_Ce,'SeriesIndex',1,'LineStyle','--','LineWidth',1.8)
plot(IRF_Cy,'SeriesIndex',1,'LineWidth',1.8,'LineStyle',':')

legend('$w_e \ne w_y$ - Expans.','$w_e \ne w_y$ - Prod.','Location','southeast','AutoUpdate','off','Interpreter','latex')

plot(IRF_C,'SeriesIndex',1,'LineWidth',1.8)
plot(IRF_YN(end-mpar_YN.oc+12,1:end-1)*100,'SeriesIndex',2,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

title('Consumption')


ylabel('PP')


subplot(4,2,4)

hold on
plot(IRF_I,'SeriesIndex',1,'LineWidth',1.8)

plot(IRF_YN(end-mpar_YN.oc+10,1:end-1)*100,'SeriesIndex',2,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

title('Investment')


subplot(4,2,5)

plot(IRF_N,'LineWidth',1.8)

hold on
plot(IRF_N_YN,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

ylabel('PP')

title('Hours')


subplot(4,2,6)

plot(IRF_W,'LineWidth',1.8)

hold on
plot(IRF_W_YN,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

title('Real wages')



subplot(4,2,7)



plot(100*(IRF_LS-par.Labshr),'LineWidth',1.8)

hold on
plot(100*(IRF_LS_YN-par.Labshr),'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';

ylabel('PP')
title('Labour share')


subplot(4,2,8)



plot(IRF_Y-IRF_N,'LineWidth',1.8)

hold on

plot(IRF_Y_YN-IRF_N_YN,'LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';


title('Labour productivity')

%%

figure(4)

subplot(1,3,1)

hold on 
%plot(Deficit,'Color',YN2Wcol,'LineWidth',1.8)
hline=refline(0,0);
hline.Color='black';

title('Debt to GDP')
ylabel('Pct of annualised GDP')

subplot(1,3,2)

hold on 
%plot(IRF_A,'Color',YN2Wcol,'LineWidth',1.8)
hline=refline(0,0);
hline.Color='black';

title('Value of fund')
ylabel('PP')

subplot(1,3,3)

hold on 
%plot(IRF_ra,'Color',YN2Wcol,'LineWidth',1.8)
hline=refline(0,0);
hline.Color='black';

title('\textbf{Illiquid return} $r_a$','Interpreter','latex')
ylabel('PP')
legend('NK','NK-YN','Location','southeast')

%legend('NK','NK-YN','NK-YN','Location','southeast')

%%
figure(7)

hold on 
%plot(IRF_ra,'Color',YN2Wcol,'LineWidth',1.8)
hline=refline(0,0);
hline.Color='black';

title('Consumption GINI coeffecient')
ylabel('')
legend('NK','NK-YN','Location','northeast')

%%

figure(3)

h = gcf;
set(h,'Units','Inches');
h.Position(3) =h.Position(3)*1.3;
h.Position(4) =h.Position(4)*1.3;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_HANK','-dpdf','-r0')

figure(4)

h = gcf;
set(h,'Units','Inches');
h.Position(3) =12;
h.Position(4) =5;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_savings','-dpdf','-r0')

figure(5)

h = gcf;
set(h,'Units','Inches');
h.Position(3) =h.Position(3)*1.3;
h.Position(4) =h.Position(4)*1.3;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_HANK2W','-dpdf','-r0')

figure(6)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
h.Position(3) =h.Position(3)*1.2;
h.Position(4) =h.Position(4)*0.9;
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_cIRFS','-dpdf','-r0')


figure(7)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
%h.Position(3) =h.Position(3)*1.2;
%h.Position(4) =h.Position(4)*0.9;
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_CGINI','-dpdf','-r0')
