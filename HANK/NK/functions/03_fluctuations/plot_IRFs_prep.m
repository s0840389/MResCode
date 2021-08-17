
IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-mpar.os,1:mpar.maxlag);

% preparation
IRF_H=100*grid.h(1:end-1)*IRF_distr(mpar.nm+mpar.nk+(1:mpar.nh-1),2:end)/par.H;
K=grid.k*IRF_distr(mpar.nm+(1:mpar.nk),:)+grid.K;
%I=(K(2:end)-(1-par.delta)*K(1:end-1));
%IRF_I=100*(I/(par.delta*grid.K)-1);

I=exp(IRF_state_sparse(end-mpar.oc+10,1:end-1))*targets.Inv;

IRF_tau=100*IRF_state_sparse(mpar.numstates-10,2:end);
IRF_K=100*IRF_state_sparse(mpar.numstates-9,1:end-1);
IRF_qs=100*IRF_state_sparse(mpar.numstates-7,2:end);
IRF_qk=100*IRF_state_sparse(mpar.numstates-8,2:end);
IRF_W=100*IRF_state_sparse(mpar.numstates-5,2:end);
IRF_G=100*IRF_state_sparse(mpar.numstates-4,1:end-1);
RB=par.RB+(IRF_state_sparse(mpar.numstates-1,2:end));
IRF_RB=100*100*(RB-par.RB);
IRF_S=100*IRF_state_sparse(end,1:end-1);

%Yss=[invmutil(mutil_c(:)); invmutil(Vk(:)); log(par.Q); log(par.PI); log(Output);...
 %   log(par.G); log(par.W)*0 ; log(par.R); log(par.PROFITS); log(par.N); log(targets.T);...
  %  log(targets.B) ;log(targets.Inv); log(par.R); log(targets.C) ];

IRF_Eq=100*IRF_state_sparse(end-mpar.oc+1,1:end-1);
IRF_PI=100*100*IRF_state_sparse(end-mpar.oc+2,1:end-1);
IRF_Y=100*IRF_state_sparse(end-mpar.oc+3,1:end-1);
IRF_G=100*IRF_state_sparse(end-mpar.oc+4,1:end-1);
IRF_PIw=100*IRF_state_sparse(end-mpar.oc+4,1:end-1);
IRF_R=100*IRF_state_sparse(end-mpar.oc+5,1:end-1);
IRF_PId=100*IRF_state_sparse(end-mpar.oc+6,1:end-1);
IRF_N=100*IRF_state_sparse(end-mpar.oc+7,1:end-1);
IRF_B=100*IRF_state_sparse(end-mpar.oc+8,1:end-1);
IRF_I=100*IRF_state_sparse(end-mpar.oc+9,1:end-1);
IRF_ra=100*IRF_state_sparse(end-mpar.oc+10,1:end-1);
IRF_C=100*IRF_state_sparse(end-mpar.oc+11,1:end-1);
IRF_giniC=100*IRF_state_sparse(end-mpar.oc+12,1:end-1);

%IRF_ly=100*IRF_state_sparse(end-mpar.oc+13,1:end-1);
%IRF_le=100*IRF_state_sparse(end-mpar.oc+14,1:end-1);
%IRF_Mg=100*IRF_state_sparse(end-mpar.oc+15,1:end-1);



IRF_M=100*grid.m*IRF_distr((1:mpar.nm),2:end)./(targets.B+par.ABS*grid.K);
M=grid.m*IRF_distr((1:mpar.nm),1:end)+targets.B-par.ABS*(K(1:end)-grid.K);

%IRF_C=100*((Y-G-I)./(Output-par.G-par.delta*grid.K)-1);


Y=Output*(1+IRF_state_sparse(end-mpar.oc+3,1:end-1));
G=par.G*(1+IRF_state_sparse(end-mpar.oc+4,1:end-1));

PI=1+IRF_state_sparse(end-mpar.oc+2,1:end-1);
Q=par.Q*(1+IRF_state_sparse(end-mpar.oc+1,1:end-1));
R=par.R*(1+IRF_state_sparse(end-mpar.oc+6,1:end-1));

IRF_RBREAL=100*100*(RB./PI-par.RB);
IRF_Q=100*100*(Q-par.Q);
IRF_D=100*100*((1+IRF_R/100)*par.R-par.R);
Deficit=100*(M(2:end) - M(1:end-1)./PI)./Y;
IRF_LP=100*100*(((Q(2:end)+R(2:end))./Q(1:end-1)-RB(1:end-1)./PI(2:end))-((1+par.R/par.Q)-par.RB));
IRF_LS=(exp(IRF_N/100)'*grid.N).*(exp(IRF_W/100)'*W_fc)./(Output*exp(IRF_Y/100)');

IRF_A=100*((qs*exp(IRF_qs/100)+par.Q*grid.K*exp(IRF_qk/100).*exp(IRF_state_sparse(mpar.numstates-6,2:end)))/targets.A - 1);

%% consumpton by wealth

IRF_hhc=IRF_state_sparse(end-length(targets.cinds)+1:end,:);

figure(5)
clf
subplot(1,3,1)
plot(100*IRF_hhc(1:6,:)')
legend(string(targets.wpcts))


subplot(1,3,2)
plot(100*IRF_hhc(7:12,:)')
legend(string(targets.wpcts))


subplot(1,3,3)
plot(100*IRF_hhc(13:18,:)')
legend(string(targets.wpcts))

%% income irf's
NN=mpar.nh*mpar.nk*mpar.nm;

tirf=length(IRF_Y);
y_irf=zeros(mpar.nm,mpar.nk,mpar.nh,tirf);


incss=par.W*par.tau/par.H*meshes.h+meshes.m(:,:,:)*(par.RB-1) +meshes.k(:,:,:)*(par.R) +meshes.m(:,:,:).*(meshes.m(:,:,:)<0)*par.borrwedge;
  
incss=incss(:);

for t=1:tirf


  WW=par.W*exp(IRF_W(t)/100)*par.tau.*exp(IRF_tau(t)/100).*exp(IRF_N(t)/100)/par.H;
  
  y_irf(:,:,:,t)= WW*meshes.h(:,:,:) +meshes.m(:,:,:)*(par.RB-1 + IRF_state_sparse(mpar.numstates-1,t))/(1+IRF_PI(:,t)/10000) +meshes.k(:,:,:)*(par.R*exp(IRF_ra(t)/100)) +meshes.m(:,:,:).*(meshes.m(:,:,:)<0)*par.borrwedge;
  

end

x=css(:);

xc=exp(IRF_hhc(1:6,:)).*x(targets.cinds(1:6));

x=reshape(y_irf,NN,tirf);

s_irf=xc./x(targets.cinds(1:6));


%% Ginis

IRF_giniw=zeros(mpar.maxlag-1,1); % wealth

IRF_top10=zeros(mpar.maxlag-1,1); % wealth top 10


for i=1:mpar.maxlag-1
%IRF_giniw(i)=100*(networthgini(IRF_distr(:,i+1)+Xss(1:mpar.nm+mpar.nk+mpar.nh),mesh,mpar,Copula)-targets.GiniW);

%IRF_giniC(i)=100*(consumptiongini(IRF_distr(:,i+1)+Xss(1:mpar.nm+mpar.nk+mpar.nh),mpar,Copula,c_irfl(:,:,:,i))-targets.GiniC);

%IRF_top10(i)=100*(top10W(IRF_distr(:,i+1)+Xss(1:mpar.nm+mpar.nk+mpar.nh),mesh,mpar,Copula)-targets.Top10);
end


% Net worth Gini
function [g]=networthgini(distr,mesh,mpar,Copula)

mplusk=mesh.k(:)+mesh.m(:);

[mplusk, IX]       = sort(mplusk);

distr2= zeros(mpar.nm+1,mpar.nk+1,mpar.nh+1);
distr2(2:end,2:end,2:end)=Copula({cumsum(distr(1:mpar.nm)),cumsum(distr(mpar.nm+1:mpar.nm+mpar.nk)),cumsum(distr(mpar.nk+mpar.nm+1:mpar.nk+mpar.nm+mpar.nh))});

distr3=diff(diff(diff(distr2,1,1),1,2),1,3);

JDredux=sum(distr3,3);

moneycapital_pdf   = JDredux(IX);

S                  = cumsum(moneycapital_pdf.*mplusk)';
S                  = [0 S];
g      = 1-(sum(moneycapital_pdf.*(S(1:end-1)+S(2:end))')/S(end));

end

% Consumption Gini

function [g]=consumptiongini(distr,mpar,Copula,c_irfl)

NN=size(c_irfl);

mplusk=squeeze(c_irfl);
mplusk=mplusk(:);

[mplusk, IX]       = sort(mplusk);

distr2= zeros(mpar.nm+1,mpar.nk+1,mpar.nh+1);
distr2(2:end,2:end,2:end)=Copula({cumsum(distr(1:mpar.nm)),cumsum(distr(mpar.nm+1:mpar.nm+mpar.nk)),cumsum(distr(mpar.nk+mpar.nm+1:mpar.nk+mpar.nm+mpar.nh))});

distr3=diff(diff(diff(distr2,1,1),1,2),1,3);

%JDredux=sum(distr3,3);

moneycapital_pdf   = distr3(IX);

S                  = cumsum(moneycapital_pdf.*mplusk)';
S                  = [0 S];
g      = 1-(sum(moneycapital_pdf.*(S(1:end-1)+S(2:end))')/S(end));

end


% top 10 percent

function [g]=top10W(distr,mesh,mpar,Copula)

mplusk=mesh.k(:)+mesh.m(:);

[mplusk, IX]       = sort(mplusk);

distr2= zeros(mpar.nm+1,mpar.nk+1,mpar.nh+1);
distr2(2:end,2:end,2:end)=Copula({cumsum(distr(1:mpar.nm)),cumsum(distr(mpar.nm+1:mpar.nm+mpar.nk)),cumsum(distr(mpar.nk+mpar.nm+1:mpar.nk+mpar.nm+mpar.nh))});

distr3=diff(diff(diff(distr2,1,1),1,2),1,3);

JDredux=sum(distr3,3);

moneycapital_pdf   = JDredux(IX);

g                  = sum(moneycapital_pdf.*mplusk.*(cumsum(moneycapital_pdf)>=0.9));
g=g/sum(moneycapital_pdf.*mplusk);


end