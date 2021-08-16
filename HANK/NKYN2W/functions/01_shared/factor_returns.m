function [ N,R_fc,Wy_fc,We_fc,Profits_fc,WW,RR,RBRB,Y,qs,par,ly,le,Mg ] = factor_returns(meshes,grid,par,mpar)
%factor_returns
  
    mc =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;
    
  Mg=1; % measure of goods

  le=mpar.MeasE; % expansionary labour
  ly=mpar.MeasY; % production labour

  N=ly+le;

  par.zn=Mg/le^par.thetan;
  
    ki=grid.K/Mg;
    lyi=ly/Mg;

  R_fc=mc*par.alphay*ki^(par.thetay*par.alphay-1)*lyi^((1-par.alphay)*par.thetay)-par.delta;
  
  Wy_fc = (1-par.alphay)*mc*ki^(par.thetay*par.alphay)*lyi^((1-par.alphay)*par.thetay-1); % wages
    
  yi = Mg*((lyi)^(1-par.alphay)*ki^(par.alphay))^par.thetay; % firm output
   
  Y=yi*Mg;

  
  We_fc=yi*(1-mc)*par.zn*par.thetan*le^(par.thetan-1);
  
  piy=Mg*(mc*yi-Wy_fc*lyi-(R_fc+par.delta)*ki); % dividends of production
  pin=Mg*yi*(1-mc)-le*We_fc;
  
  Profits_fc=piy+pin;
  
  Htoty=ly*Mg/(par.Hy*mpar.MeasY);
  Htote=le/(par.He*mpar.MeasE);
  
 NWy=Htoty*Wy_fc; % little adjustment because in egm update procedure household picks x=c-G(h,n)
  
  NWe=Htote*We_fc;
  
  WW=ones(mpar.nm,mpar.nk,mpar.nh); %Wages
  WW(:,:,1:mpar.nhy)=NWy;
  WW(:,:,mpar.nhy+1:end)=NWe;
  
  RR = R_fc; %Rental rates
  
  RBRB = par.RB/par.PI + (meshes.m<0).*(par.borrwedge/par.PI);

  qs=Profits_fc/R_fc; % stock price
  
end

