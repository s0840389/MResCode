function [ N,R_fc,W_fc,Profits_fc,WW,RR,RBRB,Y,qs,par ] = factor_returns(meshes,grid,par,mpar)
%factor_returns
  
    mc =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;
    
N=1;

  R_fc=mc*par.alphay*grid.K^(par.thetay*par.alphay-1)*N^((1-par.alphay)*par.thetay)-par.delta;
  
  W_fc = (1-par.alphay)*mc*grid.K^(par.thetay*par.alphay)*N^((1-par.alphay)*par.thetay-1); % wages
    
  Y = ((N)^(1-par.alphay)*grid.K^(par.alphay))^par.thetay; % firm output
   
  piy=(Y-W_fc*N-(R_fc+par.delta)*grid.K); % dividends of production
  pin=0;
  
  Profits_fc=piy+pin;
  
  Htot=N/par.H;
  
  par.kappaH=W_fc*par.H*par.tau/Htot^par.gamma; % work out kappa in utilty fn such that H=H
  
  NW=Htot*W_fc; % little adjustment because in egm update procedure household picks x=c-G(h,n)
  
  WW=NW*ones(mpar.nm,mpar.nk,mpar.nh);%;-par.kappaH/(par.tau*(1+par.gamma))*1./meshes.h*Htot.^(1+par.gamma); %Wages with GHH adjustment
  
 % WW(:,:,end)=Profits_fc*par.profitshare*(1-par.lumpshare);
  
  RR = R_fc; %Rental rates
  
  RBRB = par.RB/par.PI + (meshes.m<0).*(par.borrwedge/par.PI);

  qs=Profits_fc/R_fc; % stock price
  
end

