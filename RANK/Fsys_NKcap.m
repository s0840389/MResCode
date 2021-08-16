function [Difference,LHS,RHS]  = Fsys(x,xminus,y,yminus,xss,yss,p)

%xss=[sstate.K; sstate.qk; sstate.q; sstate.a; sstate.int; sstate.w; sstate.Invstate 0];

%yss=[sstate.pit; sstate.pitw; sstate.mc; sstate.N;
 %    sstate.PId; sstate.G;
 %   sstate.Inv; sstate.ra; sstate.C; sstate.pk;
 %   sstate.sy; sstate.qk; sstate.Y; sstate.d];

  yminus=exp(yss+yminus);
  y=exp(yss+y);
  xminus=exp(xss+xminus);
  x=exp(xss+x);
 

  
    controls=char('pit','pitw','mc','N','PId','G','Inv','ra','C','pk','sy','Eqk','Y','d','Ccap','Cw');
    states=char('K','qk','q','a','int','W','Invstate','taul','eint');
    delog=char('int','pk','ra','pit','d','sy','eint','pitw','taul');

% controls t


for i=1:p.numcontrols
    
    eval(sprintf('%sminus = yminus(%f);',strtrim(controls(i,:)),i))
end

% controls t+1

for i=1:p.numcontrols
    eval(sprintf('%s = y(%f);',strtrim(controls(i,:)),i))
end

% state t-1

for i=1:p.numstates
    eval(sprintf('%sminus = xminus(%f);',strtrim(states(i,:)),i))
end

% state t

for i=1:p.numstates
    eval(sprintf('%s = x(%f);',strtrim(states(i,:)),i))
end


% variabels needing logging

for i=1:size(delog,1)
        eval(sprintf('%s = log(%s);',strtrim(delog(i,:)),delog(i,:)))
        eval(sprintf('%sminus = log(%sminus);',strtrim(delog(i,:)),strtrim(delog(i,:))))
end

 % other
    wss=exp(xss(6));
    Nss=exp(yss(4));
    Gss=exp(yss(6));
    Css=exp(yss(9));
    Yss=exp(yss(13));
    Bss=exp(xss(9));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LHS=zeros(p.numcontrols+p.numstates,1);
RHS=LHS;

% (1) Euler Bonds

LHS(1)=Cwminus^-p.sigma;
RHS(1)=p.beta*(1+int)/(1+pit)*Cw^-p.sigma;

 % (2) Euler Capital

 LHS(2)=Ccapminus^(-0);
 RHS(2)=p.betaCap*(1+ra)*Ccap^(-0);
 
 % (3) Wages / Labour supply
 
 muwhat=p.psi*log(Nminus/Nss) - log(W/wss);
 
 LHS(3)=pitwminus;
 RHS(3)=p.beta*pitw+p.kappa_w*(muwhat);%+p.sigma*log(Cminus/Css);
 
 if p.Nw==0 % flexible wages
 LHS(3)=log(W/wss);
 RHS(3)=p.psi*log(Nminus/Nss);
 end
 
 
 % (4) production
 
 LHS(4)=Yminus;
 RHS(4)=(Kminus^p.alpha*(Nminus)^(1-p.alpha))^p.thetay;
 
 % (5) philips curve
 
 LHS(5)=pitminus;
 RHS(5)=p.beta*pit+p.kappa_p*(log(p.mu_p)+log(mcminus));
 
 % (6) labour demand
 
 LHS(6)=W;
 RHS(6)=mcminus*p.thetay*(1-p.alpha)*Yminus/Nminus;

 % (7) capital demand
 
 LHS(7)=pkminus;
 RHS(7)=mcminus*p.thetay*p.alpha*Yminus/Kminus;
 
%  (8) Dividend
 
 divy=Yminus-W*Nminus-Kminus*pkminus;
 
 LHS(8)=PIdminus;
 RHS(8)=divy;
 
% (9) capital price

 LHS(9)=qk;
% RHS(9)=1+p.tau*(Invminus/Kminus-p.delta);
 RHS(9)=1+p.tau*log(Invminus/Invstateminus)+p.tau/2*log(Invminus/Invstateminus)^2- p.tau/(1+ra)*Inv/Invminus*log(Inv/Invminus);
% (10) expected capital price
 
 LHS(10)=Eqkminus;
 RHS(10)=qk;
 
% (11) expected capital price
 
 LHS(11)=Eqkminus;
 %RHS(11)=(1/(1+ra))*(pk-p.tau/2*(Inv/K-p.delta)^2+p.tau*(Inv/K-p.delta)*Inv/K+(1-p.delta)*Eqk);
  RHS(11)=1/(1+ra)*(pk+(1-p.delta)*Eqk);
 % (12) capital lom
 
 LHS(12)=K;
% RHS(12)=(1-p.delta)*Kminus+Invminus-p.tau/2*(Invminus/Kminus-p.delta)^2*Kminus;
  RHS(12)=(1-p.delta)*Kminus+Invminus-p.tau/2*log(Invminus/Invstateminus)^2*Invminus;

% (13) fund value

 LHS(13)=a;
 RHS(13)=K*qk+q;

% (14) ra
 
 LHS(14)=raminus*aminus;
 RHS(14)=Kminus*qk*(1-p.delta)-Kminus*qkminus+pkminus*Kminus+PIdminus +q-qminus;

% (15) fund lom
 
 LHS(15)=a;
 RHS(15)=(1+raminus)*aminus+dminus;
 
% (16) Government spending

LHS(16)=Gminus;
RHS(16)=Gss;
 

% (17) resource constraint
 
 LHS(17)=Yminus;
 RHS(17)=Cminus + Gminus + Invminus + p.tau/2*Invminus*log(Invminus/Invstateminus)^2;
 
% (18) Interest rate

 LHS(18)=int;
 RHS(18)=p.rhoint*intminus+ (1-p.rhoint)*(p.phipi*(pitminus-p.pitstar) - p.phiy*log(Yminus/Yss)+ p.intstar) + eintminus;
 
% (19) wage inflation


    LHS(19)=pitw;
    RHS(19)=log(W/Wminus)+pitminus;

    if p.Nw==0
       LHS(19)=pitwminus;
       RHS(19)=log(W/Wminus);
    end
    
 % (20) labor share
 
 LHS(20)=syminus;
 RHS(20)=W*Nminus/Yminus;

 % (21) Investment state variable
 LHS(21)=Invstate;
 RHS(21)=Invminus;
 

 % (22) Debt

LHS(22)       = (Gminus);
RHS(22) =   W*Nminus*(1-taul);


% (23) total consumption
 
 LHS(23)=Cminus;
 RHS(23)=p.shrCap*Ccapminus+(1-p.shrCap)*Cwminus;


 % (24) capatalsit consumption
 
 LHS(24)=Ccapminus;
 RHS(24)=-dminus/p.shrCap;
 
 % (25 shock)
 
 RHS(25)=eint;
 LHS(25)=0*eintminus;

%% Difference
Difference=((LHS-RHS));

end
