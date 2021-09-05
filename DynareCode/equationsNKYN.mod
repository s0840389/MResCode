// Specifies equations for firstOrderDynamics.mod


# CTREND=1+ctrend;

// Euler bonds (1)

	exp(c)^-ssigma=bbeta*CTREND^(-ssigma)*(1+int+intstar)/(1+pitstar+pit(+1))*exp(c(+1))^-ssigma;

// Euler capital (2)

	exp(c)^-ssigma=bbeta*CTREND^(-ssigma)*(1+ra(+1))*exp(c(+1))^-ssigma;

// Labour supply (3)

 # muwhat=ppsi*(l-lss) - (w-wss);
 
 pitw= bbeta/CTREND^ssigma*pitw(+1) +kappa_w*(muwhat)+zwmk;
 
// Production y (4)

	exp(yi)=exp(zy)*(exp(k(-1))/exp(N))^aalpha*exp(ly)^(1-aalpha);

// Production N (4b)

	exp(N)=exp(zy)*zn*exp(le)^tthetan;

// Philips Curve (5)
 
	pit=bbeta/CTREND^ssigma*pit(+1) +kappa_p*(log(mu_p)+mc)+zmk;

 // Labour demand production (6)
 
	exp(w)=exp(mc)*(1-aalpha)*exp(yi)/exp(ly);
 
 // capital demand (7)
 
	pk=exp(mc)*aalpha*exp(yi)/(exp(k(-1))/exp(N));
 
 // Labour demand expansion (10)
 
 # piN=exp(yi)*(1-exp(mc));
 
	exp(w)=piN*exp(zy)*zn*tthetan*exp(le)^(tthetan-1);
 
	exp(l)=exp(N)*exp(ly)+exp(le);
 
	exp(y)=exp(yi)*exp(N);
 
// dividend (11)
 
 # divy = exp(N)* (exp(yi)*exp(mc)	-	exp(w)*exp(ly)	- exp(k(-1))/exp(N)*pk);
 # divn = exp(N) *piN - exp(w)*exp(le);
 
// dividend (8)
 
	exp(PId)=	divy+divn;
	
// capital price (9)
 
 exp(qk)*exp(zi)=1+ttau*log(exp(Inv)/exp(Inv(-1))) + ttau/2*log(exp(Inv)/exp(Inv(-1)))^2 - ttau/(1+ra(+1))*exp(Inv(+1))/exp(Inv)*log(exp(Inv(+1))/exp(Inv));
 
// capital return (10)

	exp(qk)=1/(1+ra(+1))*(pk(+1)+(1-ddelta)*exp(qk(+1))-zrp);
 
// capital lom (11)
 
	CTREND*exp(k)=(1-ddelta)*exp(k(-1))+exp(Inv)*exp(zi) -ttau/2*log(exp(Inv)/exp(Inv(-1)))^2*exp(Inv);
 
// fund value (12)
 
	exp(a)=exp(k)*exp(qk)+exp(q);
 
// fund lom (13-14) 
  
	exp(a(-1))*ra=exp(k(-1))*exp(qk)*(1-ddelta)-exp(k(-1))*exp(qk(-1))+ pk*exp(k(-1))+exp(PId) + exp(q)-exp(q(-1));
  
	CTREND*exp(a)=exp(a(-1))*(1+ra)+d;
 
// government spending (15)
 
	g=(1-rhogov)*gss+rhogov*g(-1)+egov;

// resource constraint (16)

	exp(y)=exp(c)+exp(g)+exp(Inv)+ttau/2*exp(Inv)*log(exp(Inv)/exp(Inv(-1)))^2;
	
// interest rate (17)
 
	int=rhoint*int(-1)+ (1-rhoint)*(phipi*(pit)+phiy*(y-yf)+phiyy*(y-yf-y(-1)+yf(-1))) + eint;
 
// DEBT (18)

# taxrev=exp(w)*exp(l)*(1-taxrate);

CTREND*exp(B)=(1+int(-1))/(1+pit)*exp(B(-1))+exp(g)-taxrev;
 
// taxes (19)

(1-taxrate)/(1-taxss)= ((1-taxrate(-1))/(1-taxss))^rho_tax *((exp(B)/exp(Bss))^gamma_taxB*(exp(y)/exp(yss))^gamma_taxY)^(1-rho_tax) ;

// wage inflation (19)

	pitw=pit+w-w(-1);

// labour share (20)

	//sy=exp(w)*exp(l)/(exp(y));
	  
// shocks
	
	zy=rhozy*zy(-1)+ezy; // tfp
	zi=rhozi*zi(-1)+ezi; // ist
	zmk=rhozmk*zmk(-1)+ezmk; // mark up
	zwmk=rhozwmk*zwmk(-1)+ezwmk; // wage mark up
	zrp=rhozrp*zrp(-1)+ezrp; // risk premium
	 
	 
// flexible economy [ yf cf lyf wf pkf qkf kf Invf yif Nf lef lf]

// Euler capital (1f)

	exp(qk)*exp(cf)^-ssigma=bbeta*CTREND^(-ssigma)*((1-ddelta)*exp(qk(+1))+pk(+1))*exp(c(+1))^-ssigma;

// Labour supply (2f)

 wf-wss=ppsi*(lf-lss);
 
// Production y (3f)

	exp(yif)=exp(zy)*(exp(kf(-1))/exp(Nf))^aalpha*exp(lyf)^(1-aalpha);

 // Labour demand production (4f)
 
	exp(wf)=exp(mcss)*(1-aalpha)*exp(yif)/exp(lyf);
 
 // capital demand (5f)
 
	pkf=exp(mcss)*aalpha*exp(yif)/(exp(kf(-1))/exp(Nf));
  
// capital price (6f)
 
 exp(qkf)*exp(zi)=1+ttau*log(exp(Invf)/exp(Invf(-1))) + ttau/2*log(exp(Invf)/exp(Invf(-1)))^2 - bbeta/CTREND^ssigma*ttau/(1+ra(+1))*exp(Inv(+1))/exp(Inv)*log(exp(Inv(+1))/exp(Inv));
 
// capital lom (7f)
 
	CTREND*exp(kf)=(1-ddelta)*exp(kf(-1))+exp(Inv)*exp(zi) -ttau/2*log(exp(Invf)/exp(Invf(-1)))^2*exp(Invf);
 
// resource constraint (8f)

	exp(yf)=exp(cf)+exp(gss)+exp(Invf)+ttau/2*exp(Invf)*log(exp(Invf)/exp(Invf(-1)))^2;
	
	exp(yf)=exp(Nf)*exp(yif);
	
 #piNf=exp(yif)*(1-exp(mcss));
 
	exp(wf)=piNf*exp(zy)*zn*tthetan*exp(lef)^(tthetan-1);
 
	exp(lf)=exp(Nf)*exp(lyf)+exp(lef);
	
	exp(Nf)=exp(zy)*zn*exp(lef)^tthetan;

	
	 
// measurment equations

	dy=100*(y-y(-1))+100*ctrend;
	dc=100*(c-c(-1))+100*ctrend;
	dinve=100*(Inv-Inv(-1))+100*ctrend;
	dw=100*(w-w(-1))+100*ctrend;
	pinfobs = 100*(pit)+100*pitstar;
	robs =    100*(int)+100*intstar;
	labobs = (ly-lss)*100 + constelab;
	
	
	sy=100*0.62*(w-wss+l-lss-y+yss) +constelab;
	