

load('../../steadystates/YNfund_60_15.mat')



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
    x=hx*x;
end

x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse_rb=[];
x=x0;

% y=ybar+gx*(x-xbar)
% x'=xbar+hx(x-xbar)

hx_rb=hx;

hx_rb(1:end-2,:)=0;

hx_rb(end-1,1:end-2)=0;

for t=1:mpar.maxlag
    IRF_state_sparse_rb(:,t)=(MX*x)';
    x=hx_rb*x;
end

figure(1)
clf
plot(IRF_state_sparse(end-mpar.oc+12,:))
hold on
plot(IRF_state_sparse_rb(end-mpar.oc+12,:),'r--')
