function [Difference,LHS,RHS]  = Fsys(x,xminus,y,yminus,xss,yss,p)

%xss=[sstate.K; sstate.qk; sstate.q; sstate.a; sstate.int; sstate.w; 0];


%yss=[sstate.pit; sstate.pitw; sstate.mc; sstate.N;
 %    sstate.PId; sstate.G; 
  %  sstate.Inv; sstate.ra; sstate.C; sstate.pk;
  %  sstate.sy; sstate.qk; sstate.Y; sstate.d;
  %  sstate.Mg; sstate.le; sstate.ly; sstate.ki; sstate.yi];

  yminus=exp(yss+yminus);
  y=exp(yss+y);
  xminus=exp(xss+xminus);
  x=exp(xss+x);
 
    controls=char('pit','pitw','mc','N','PId','G','Inv','ra','C','pk','sy','Eqk','Y','d','Mg','le','ly','ki','yi','Ccap','Cw');
    states=char('K','qk','q','a','int','W','taul','Invstate','eint');
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
    Bss=exp(xss(8));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LHS=zeros(p.numcontrols+p.numstates,1);
RHS=LHS;

% (1) Euler Bonds

LHS(1)=Cwminus^-p.sigma;
RHS(1)=p.beta*(1+int)/(1+pit)*Cw^-p.sigma;

 % (2) Euler Capital
 

 LHS(2)=Ccap^0;
 RHS(2)=p.betaCap*(1+ra)*C^-0;
 
 % (3) Wages / Labour supply
 
 muwhat=p.psi*log(Nminus/Nss) - log(W/wss);
 
 LHS(3)=pitwminus;
 RHS(3)=p.beta*pitw+p.kappa_w*(muwhat);%+p.sigma*log(Cminus/Css);
 
 if p.Nw==0 % flexible wages
 LHS(3)=log(W/wss);
 RHS(3)=p.psi*log(Nminus/Nss);
 end
 
 
 % (4) production
 
 LHS(4)=yiminus;
 RHS(4)=(kiminus^p.alphay*(lyminus)^(1-p.alphay))^p.thetay;
  
 % (5) Measure of goods
 
 LHS(5)=Mgminus;
 RHS(5)=p.zn*(leminus)^p.thetan;
 
 % (6) philips curve
 
 LHS(6)=pitminus;
 RHS(6)=p.beta*pit+p.kappa_p*(log(p.mu_p)+log(mcminus));
 
 % (7) production labour demand
 
 LHS(7)=W;
 RHS(7)=mcminus*p.thetay*(1-p.alphay)*yiminus/lyminus;

  % (8) expansion labour demand
 
 LHS(8)=W;
 RHS(8)=yiminus*(1-mcminus)*p.thetan*Mgminus/leminus;

 % (9) capital demand
 
 LHS(9)=pkminus;
 RHS(9)=mcminus*p.thetay*p.alphay*yiminus/kiminus;
 
%  (10) Dividend
 
 divy=Mgminus*(yiminus*mcminus-W*lyminus-kiminus*pkminus);
 divn=Mgminus*yiminus*(1-mcminus)-W*leminus;
 
 LHS(10)=PIdminus;
 RHS(10)=divy+divn;
 
% (11) capital price

  LHS(11)=qk;
  RHS(11)=1+p.tau*log(Invminus/Invstateminus)+p.tau/2*log(Invminus/Invstateminus)^2- p.tau/(1+ra)*Inv/Invminus*log(Inv/Invminus);

% (12) expected capital price
 
 LHS(12)=Eqkminus;
 RHS(12)=qk;
 
% (13) expected capital price
 
 LHS(13)=Eqkminus;
 RHS(13)=1/(1+ra)*(pk+(1-p.delta)*Eqk);
  
 % (14) capital lom
 
 LHS(14)=K;
 RHS(14)=(1-p.delta)*Kminus+Invminus-p.tau/2*log(Invminus/Invstateminus)^2*Invminus;
 
% (15) fund value

 LHS(15)=a;
 RHS(15)=K*qk+q;

% (14) ra
 
 LHS(16)=raminus*aminus;
 RHS(16)=Kminus*qk*(1-p.delta)-Kminus*qkminus+pkminus*Kminus+PIdminus +q-qminus;

% (17) fund lom
 
 LHS(17)=a;
 RHS(17)=(1+raminus)*aminus+dminus;
 
% (18) Government spending

LHS(18)=Gminus;
RHS(18)=(1-taul)*Nminus*W;
 
% (19) resource constraint
 
 LHS(19)=Yminus;
 RHS(19)=Cminus + Gminus + Invminus + p.tau/2*log(Invminus/Invstateminus)^2*Invminus;
 
% (20) Interest rate

 LHS(20)=int;
 RHS(20)=p.rhoint*intminus+ (1-p.rhoint)*(p.phipi*(pitminus-p.pitstar) - p.phiy*log(Yminus/Yss)+ p.intstar) + eintminus;
 
% (21) wage inflation


    LHS(21)=pitw;
    RHS(21)=log(W/Wminus)+pitminus;

    if p.Nw==0
       LHS(19)=pitwminus;
       RHS(19)=log(W/Wminus);
    end
    
 % (22) labor share
 
 LHS(22)=syminus;
 RHS(22)=W*Nminus/Yminus;

 % (23) Capital

 LHS(23)=Kminus;
 RHS(23)=Mgminus*kiminus;
 
 % (24) Output
 
 LHS(24)=Yminus;
 RHS(24)=Mgminus*yiminus;
 
 % (25) Labour
 
 LHS(25)=Nminus;
 RHS(25)=leminus+lyminus*Mgminus;


% (27) taxes

    LHS(26)       = Gss;

    RHS(26) = Gminus;

    % (23) total consumption
 
 LHS(27)=Cminus;
 RHS(27)=p.shrCap*Ccapminus+(1-p.shrCap)*Cwminus;


 % (24) capatalsit consumption
 
 LHS(28)=Ccapminus;
 RHS(28)=-dminus/p.shrCap;
 
% (28) Investment state

    LHS(29)=Invstate;
    RHS(29)=Invminus;

 % (29) shock
 
 RHS(30)=eint;
 LHS(30)=0*eintminus;

%% Difference
Difference=((LHS-RHS));

end
