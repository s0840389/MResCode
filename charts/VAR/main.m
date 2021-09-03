clear all; close all; 

% Create DATASET
%%%%%%%%%%%%%%%%
DATASET.TSERIES = xlsread('VARDATA.xlsx','Quarterly2');
DATASET.LABEL   = {'DATES','R','Y','P','CPI','PoC','W','LS','Rinst'};
DATASET.VALUE   = [  1,       2,     3,     4,    5 ,    6 ,     7, 8, 9 ];
DATASET.UNIT    = [  0,       2,     2,     2,    2 ,    2,      2 ,2 ,2];


DATASET.MAP = containers.Map(DATASET.LABEL,DATASET.VALUE);

% VAR specification
%%%%%%%%%%%%%%%%%%%%
VAR.p                = 3;  % Number of Lags
VAR.irhor            = 20; % Impulse Response Horizon
VAR.select_vars      = {'R','Y','P','CPI','PoC','W','LS'};
VAR.vars             = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
VAR.MAP              = containers.Map(VAR.select_vars,1:size(VAR.vars,2));
VAR.proxies          = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,{'Rinst'})));
VAR.DET              = ones(length(VAR.vars),1); % Deterministic Terms

VAR        = doProxySVAR(VAR);

% Inference:
%%%%%%%%%%%%
% method  1: Mertens and Ravn (2013) wild bootstrap 
%         2: Montiel-Olea Stock Watson (2016) parametric bootstrap
%         3: Delta Method 
%         4: Montiel-Olea Stock Watson (2016) asy weak IV
%         5: Jentsch and Lunsford Moving Block Bootstrap
%         6: Jentsch and Lunsford Moving Block Bootstrap (adjusted to allow non zero-mean proxies)

nboot     = 5000;        % Number of Bootstrap Samples (Paper does 5000)
clevel    = 68;          % Bootstrap Percentile Shown
BlockSize = 19;          % size of blocks in the MBB bootstrap
seed      = 2;           % seed for random number generator
%rng(seed);               % iniate the random number generator

VARci_wildbs   = doProxySVARci(VAR,clevel,1,nboot); 
%VARci_mswbs    = doProxySVARci(VAR,clevel,2,nboot);
%VARci_delta    = doProxySVARci(VAR,clevel,3);
VARci_mbb      = doProxySVARci(VAR,clevel,5,nboot,BlockSize);


%%


meanLS=mean(exp(VAR.vars(:,7)/100));

figure(1)
clf

subplot(1,2,1)

lo=VARci_mbb.irsL(:,7)*-1*meanLS;
hi=VARci_mbb.irsH(:,7)*-1*meanLS;;

   h = area([1:20],[lo,hi - lo],'FaceAlpha',0.20);
            
            set(h(2),'FaceColor',[0 0 .9])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
hold on
            plot([1:20],-1*VAR.irs(:,7)*meanLS,'blue','LineWidth',1.8)


hline=refline(0,0);
hline.Color='black';
xlabel('Qtr')
ylabel('Percentage points')
title('Labour share')

xlim([1,20])


subplot(1,2,2)

lo=VARci_mbb.irsL(:,7)+VARci_mbb.irsL(:,6);
hi=VARci_mbb.irsH(:,7)+VARci_mbb.irsH(:,6);

   h = area([1:20],[lo,hi - lo],'FaceAlpha',0.20);
            
            set(h(2),'FaceColor',[0 0 .9])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
hold on
            plot([1:20],VAR.irs(:,7)+VAR.irs(:,6),'blue','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
xlabel('Qtr')
ylabel('Percent')

title('Labour productivity')
xlim([1,20])



%% appendix chart

labs={'R','log(GDP)','log(GDP.Def)','log(CPI)','log(COM)','log(Wage)','Labour share'}

figure(2)
clf

for i=1:7
    
subplot(3,3,i)

hold on


if i==1
lo=-1*VARci_mbb.irsL(:,i);
hi=-1*VARci_mbb.irsH(:,i);
ols=-1*VAR.irs(:,i);
end
    
if i>1 & i<7
lo=-1*VARci_mbb.irsL(:,i)*100;
hi=-1*VARci_mbb.irsH(:,i)*100;
ols=-1*VAR.irs(:,i)*100;
end


if i==7
lo=-1*VARci_mbb.irsL(:,i)*meanLS;
hi=-1*VARci_mbb.irsH(:,i)*meanLS;
ols=-1*VAR.irs(:,i)*meanLS;
end

h = area([1:20],[lo,hi - lo],'FaceAlpha',0.20);
            
            set(h(2),'FaceColor',[0 0 .9])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
hold on
            plot([1:20],ols,'blue','LineWidth',1.8)


    ylabel('Difference')
   
    
    
    title(labs(i))
end

%% saving

figure(1)
h = gcf;
set(h,'Units','Inches');
h.Position(3) =h.Position(3)*1.2;
h.Position(4) =h.Position(4)*0.8;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'../IRF_US_LS','-dpdf','-r0')


figure(2)
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'../IRF_VAR','-dpdf','-r0')


%% First stage F test

for i=1:1
X=[VAR.m];
Y=VAR.res(:,i);
lmF=fitlm(X,Y,'RobustOpts','on');
ftab=table2array(lmF.anova);
F(i)=ftab(1,4);
bb(i)=table2array(lmF.Coefficients(2,1));
end






