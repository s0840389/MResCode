clc
clear

addpath(genpath('../HANKYN/functions'))

casename            = '../../steadystates/YNfund2W_60_15.mat';

load(casename)

T=101;
Tf=24;

%% Produce IRFs
disp('Calculating IRFs.');

mpar.maxlag=T; % Quarters

x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_decomp=[];
x=x0;

% y=ybar+gx*(x-xbar)
% x'=xbar+hx(x-xbar)


for t=1:mpar.maxlag
    IRF_decomp(:,t)=(MX*x)';
    x=hx*x*exp(max(t-25,0)*-1/250);;
end


%% solve consumption

T=T-1;

Vk_decomp=zeros(mpar.nm,mpar.nk,mpar.nh,T);
mutil_c_decomp=zeros(mpar.nm,mpar.nk,mpar.nh,T);
c_decomp=zeros(mpar.nm,mpar.nk,mpar.nh,T);
m_a_star_decomp=zeros(mpar.nm,mpar.nk,mpar.nh,T);
m_n_star_decomp=zeros(mpar.nm,mpar.nk,mpar.nh,T);
k_a_star_decomp=zeros(mpar.nm,mpar.nk,mpar.nh,T);


Vk_decomp(:,:,:,end)=Vk;
mutil_c_decomp(:,:,:,end)=mutil_c;
c_decomp(:,:,:,end)=css;

% eq prices

RB=par.RB*exp(IRF_decomp(mpar.numstates-1,1:end-1));

PI=par.PI*exp(IRF_decomp(end-mpar.oc+2,1:end-1));

Wy=par.Wy*exp(IRF_decomp(mpar.numstates-3,2:end));
We=par.We*exp(IRF_decomp(mpar.numstates-2,2:end));

le=grid.le*exp(IRF_decomp(end-mpar.oc+14,1:end-1));
ly=grid.ly*exp(IRF_decomp(end-mpar.oc+13,1:end-1));


tau=par.tau*exp(IRF_decomp(mpar.numstates-8,2:end));

ra=par.R*exp(IRF_decomp(end-mpar.oc+11,1:end-1));
Mg=grid.Mg*exp(IRF_decomp(end-mpar.oc+15,1:end-1));


f = waitbar(0,'1','Name','Computing value fns',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for tb=1:T-1

    t=T-tb;
    
  waitbar(tb/(T-1),f,sprintf('%i ',round(100*tb/(T-1),0)))
  
  Htoty=ly(t)/(par.Hy*mpar.MeasY);
  Htote=le(t)/(par.He*mpar.MeasE);
  
  NWy=Htoty*Wy(t); % little adjustment because in egm update procedure household picks x=c-G(h,n)
  NWe=Htote*We(t);
  
  WW=ones(mpar.nm,mpar.nk,mpar.nh); %Wages
  WW(:,:,1:mpar.nhy)=NWy;
  WW(:,:,mpar.nhy+1:end)=NWe;
    
% Incomes (grids)
inc.labor   = tau(t)*WW(t).*(meshes.h);%+Profitminus*par.lumpshare*par.tau;
inc.rent    = meshes.k*ra(t);%+(1-par.lumpshare)*meshes.k/Kminus*par.tau*Profitminus;
inc.capital = meshes.k;
inc.money   = (RB(t)/PI(t)).*meshes.m...
    + (meshes.m<0).*(par.borrwedge/PI(t)).*meshes.m;

% value functions
EVk = reshape(reshape(Vk_decomp(:,:,:,t+1),[mpar.nm*mpar.nk mpar.nh])*P_H',[mpar.nm,mpar.nk mpar.nh]);
RBaux = RB(t)/PI(t+1) + (meshes.m<0).*(par.borrwedge/PI(t+1));
mutilc_temp=mutil_c_decomp(:,:,:,t+1);
mutilc_temp=mutilc_temp(:);
EVm = reshape(reshape(RBaux(:).*mutilc_temp,[mpar.nm*mpar.nk mpar.nh])*P_H',[mpar.nm,mpar.nk mpar.nh]);

% solve c
[c_a_star,m_a_star,k_a_star,c_n_star,m_n_star] = EGM_policyupdate(EVm,EVk,1,PI(t),RB(t),inc,meshes,grid,par,mpar);

m_a_star_decomp(:,:,:,t)=m_a_star;
m_n_star_decomp(:,:,:,t)=m_n_star;
k_a_star_decomp(:,:,:,t)=k_a_star;

c_decomp(:,:,:,t)=par.nu*c_a_star+(1-par.nu)*c_n_star;

mutil_c_n = 1./(c_n_star.^par.xi); % marginal utility at consumption policy no adjustment
mutil_c_a = 1./(c_a_star.^par.xi); % marginal utility at consumption policy adjustment
mutil_c_decomp(:,:,:,t) = par.nu.*mutil_c_a + (1-par.nu).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)

%% Update marginal Value of Capital

EVk     = reshape(Vk_decomp(:,:,:,t+1),[mpar.nm*mpar.nk, mpar.nh])*P_H';
EVkl=reshape(EVk,[mpar.nm mpar.nk mpar.nh]);

Vk_next=zeros(mpar.nm,mpar.nk,mpar.nh);


for i=1:mpar.nh
    
Vk_nextInt = griddedInterpolant(meshes.m(:,:,i),meshes.k(:,:,i),EVkl(:,:,i));

Vk_next(:,:,i)=Vk_nextInt(m_n_star(:,:,i),meshes.k(:,:,i));

end

Vk_decomp(:,:,:,t)  = par.nu.*(ra(t)+1).*mutil_c_a + (1-par.nu).*ra(t).*mutil_c_n+ par.beta.*(1-par.nu).*Vk_next; % Expected marginal utility at consumption policy (w &w/o adjustment)


end

delete(f)

%% Forward IRF 

JD=zeros(mpar.nm,mpar.nk,mpar.nh,40);

JD(:,:,:,1)=joint_distr;


f = waitbar(0,'1','Name','Computing forward IRF',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for t=1:Tf-1
    
  waitbar(t/(Tf-1),f,sprintf('%f ',round(100*t/(Tf-1),0)))
  
    % Initialize matrices
weight11  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight21  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight22  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);

% find next smallest on-grid value for money and capital choices
[Dist_m_a,idm_a] = genweight(m_a_star_decomp(:,:,:,t),grid.m);
[Dist_m_n,idm_n] = genweight(m_n_star_decomp(:,:,:,t),grid.m);
[Dist_k,idk_a]   = genweight(k_a_star_decomp(:,:,:,t),grid.k);
idk_n = repmat(ones(mpar.nm,1)*(1:mpar.nk),[1 1 mpar.nh]); %This is the actual point on grid


% Transition matrix adjustment case
idm_a=repmat(idm_a(:),[1 mpar.nh]);
idk_a=repmat(idk_a(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nk*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:),idk_a(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:),idk_a(:)+1,idh(:));
index21 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:)+1,idk_a(:),idh(:));
index22 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_a(:)+1,idk_a(:)+1,idh(:));


for hh=1:mpar.nh
    
    %Corresponding weights
    weight21_aux =  Dist_m_a(:,:,hh).*(1-Dist_k(:,:,hh));
    weight11_aux = (1-Dist_m_a(:,:,hh)).*(1-Dist_k(:,:,hh));
    weight22_aux =  (Dist_m_a(:,:,hh)).*(Dist_k(:,:,hh));
    weight12_aux =  (1-Dist_m_a(:,:,hh)).*(Dist_k(:,:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P_H(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P_H(hh,:);
    weight21(:,:,hh)=weight21_aux(:)*P_H(hh,:);
    weight22(:,:,hh)=weight22_aux(:)*P_H(hh,:);
end
%Dimensions (mxk,h,h')
weight11=permute(weight11,[1 3 2]);
weight22=permute(weight22,[1 3 2]);
weight12=permute(weight12,[1 3 2]);
weight21=permute(weight21,[1 3 2]);
rowindex=repmat(1:mpar.nm*mpar.nk*mpar.nh,[1 4*mpar.nh]);

H_a=sparse(rowindex(:),[index11(:); index21(:); index12(:); index22(:)],...
    [weight11(:); weight21(:); weight12(:); weight22(:)],mpar.nm*mpar.nk*mpar.nh,mpar.nm*mpar.nk*mpar.nh); % mu'(h',k'), a without interest

% Policy Transition Matrix for no-adjustment case
weight11  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
weight21  = zeros(mpar.nm*mpar.nk, mpar.nh,mpar.nh);
idm_n=repmat(idm_n(:),[1 mpar.nh]);
idk_n=repmat(idk_n(:),[1 mpar.nh]);
index11 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_n(:),idk_n(:),idh(:));
index21 = sub2ind([mpar.nm mpar.nk mpar.nh],idm_n(:)+1,idk_n(:),idh(:));


for hh=1:mpar.nh
    
    %Corresponding weights
    weight21_aux =  Dist_m_n(:,:,hh);
    weight11_aux = (1-Dist_m_n(:,:,hh));
    
    weight21(:,:,hh)=weight21_aux(:)*P_H(hh,:);
    weight11(:,:,hh)=weight11_aux(:)*P_H(hh,:);
    
end
weight11=permute(weight11,[1 3 2]);
weight21=permute(weight21,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nk*mpar.nh,[1 2*mpar.nh]);

H_n=sparse(rowindex,[index11(:); index21(:)],...
    [weight11(:); weight21(:)],mpar.nm*mpar.nk*mpar.nh,mpar.nm*mpar.nk*mpar.nh); % mu'(h',k'), a without interest

% Joint transition matrix and transitions
H=par.nu*H_a + (1-par.nu)*H_n;

JDminus=JD(:,:,:,t);

JD_new=JDminus(:)'*H;
JD(:,:,:,t+1)= reshape(JD_new(:),[mpar.nm,mpar.nk,mpar.nh]);

    
end

delete(f)


IRF_Cdecomp=(sum(reshape(c_decomp(:,:,:,1:Tf),NN,Tf).*reshape(JD(:,:,:,1:Tf),NN,Tf))/targets.C-1)*100;

figure(110)
clf
plot(IRF_Cdecomp)



function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function


