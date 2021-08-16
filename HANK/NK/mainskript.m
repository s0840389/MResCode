%%======================================================
% Jamie Lenney Dissertation (Two sector HANK)

%% Initialize workspace and load directories
clc
clear

addpath(genpath('functions'))


UCL=0;

%% Select options
% Search for steady state and name it
FindNewSS           = false;
casename            = '../../../../steadystates/NKfund_60_15.mat';
filename=casename;

if UCL==1
casename            = '../steadystates/NKfund_60_15.mat';
filename=casename;
end

%% Solve for Steady state
if FindNewSS  
    
    disp('Solving Steady State by EGM')

    % Set parameters
    defineSS_pars
    
    mainskript_steadystate

    disp('Steady State solved')
else
    disp('Loading Steady State')

    load(casename)
   
end

%%
doirf=true;
mpar.overrideEigen  = true; % Warning appears, but critical Eigenvalue shifted
print_IRFs          = true; % Plot and print as PDF

if doirf

%% Select aggregate shock
% aggrshock           = 'TFP';
%aggrshock           = 'Uncertainty';
aggrshock           = 'MP';

par.rhoS            = 0.00;    % Persistence 
par.sigmaS          = 0.01;    % STD s

disp(strcat('Calculating IRF: ',aggrshock))

%% Produce matrices to reduce state-space
mainskript_statereduc

%% Initialize System
disp('Computing system for SS.');
F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,indexMUdct,indexVKdct,DC1,DC2,DC3,...
    par,mpar,grid,targets,Copula,P_H,aggrshock);

[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);

if max(abs(Fss))>0.001
    [Fss LHS RHS]
    error('Fss NE 0')
end
%% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,mpar,par);

%% Produce IRFs
disp('Calculating IRFs.');

mpar.maxlag=40; % Quarters

x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

% y=ybar+gx*(x-xbar)
% x'=xbar+hx(x-xbar)


for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x%*exp(-max(0,t-25)/250);
end

plot_IRFs_prep

if print_IRFs


irflist=char('RB','PI','Y','C','I','M','LS','N','W');

figure(100)
clf

for i=1:9
    
subplot(3,3,i)

eval(sprintf('plot(IRF_%s)',irflist(i,:)))
eval(sprintf('title("%s")',irflist(i,:)))

if irflist(i,:)=='LS'
hold on
plot(ones(size(IRF_LS))*targets.LS)
hold off
end

end

saveas(gcf,'IRF6.jpg')

figure(101)
subplot(3,1,1)
bar(moneydist)
subplot(3,1,2)
bar(capdist)
subplot(3,1,3)
plot(squeeze(c_a_guess(:,25,:)))
legend(num2str([1:mpar.nh]'))


end

end



%% Saving again (with IRFs)

save(filename)