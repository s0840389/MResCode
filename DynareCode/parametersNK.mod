// Declare parameters and load in their values for firstOrderDynamics_polynomials.mod
//

//----------------------------------------------------------------
// Preliminaries
//----------------------------------------------------------------

// Load in files containing parameters
//load('p');


//---------------------------------------------------------
// Economic parameters
//----------------------------------------------------------

parameters ssigma eeps ppsi aalpha ddelta tthetay BY GY
phipi pitstar intstar    cchi0 cchi1 cchi2 lss wss ttau  
ctrend  constelab kappa_p kappa_w taxss gss mu_p Bss yss rho_tax gamma_taxB gamma_taxY phiy
rhozy rhozi rhozmk rhozwmk rhozrp rhogov rhoint  bbeta cbeta mcss phiyy;

clab=0;

ctrend=(1.004-1);
//intstar=0.05/4;
pitstar=0.03/4; % inflation target

ssigma=1.0;	% crra

cbeta=0.0025;

bbeta=1/(1+cbeta);

intstar=(1+ctrend)^ssigma*(1+pitstar)/(bbeta)-1; % r*



eeps=20;	% epsilon
ppsi=2;	% labour elast (inverse frish)

mu_p=eeps/(eeps-1);

aalpha=1-0.62*mu_p;	% capital elasticity
ddelta=0.07/4;
tthetay=1.00;	% scale Y

mu_p=eeps/(eeps-1);

ttau=1.8;

% tax rules
rho_tax=0.55; %autocorrelation
gamma_taxB=0.777; %debt
gamma_taxY=2.646; %output

BY=1.6;
GY=0.175;

taxss=1-(GY+BY*intstar/(1+pitstar))/0.62;

phipi=1.8; % taylor rule
phiy=0.1; % taylor rule
phiyy=0.05;

kappa_p=(1-0.75)*(1-0.75*bbeta)/0.75;
kappa_w=(1-0.75)*(1-0.75*bbeta)/0.75;



% adjustment costs 
cchi0=0*0.04;
cchi1=0*0.90;
cchi2=1.3;

% trends


conster=(intstar)*100;
constelab=0;

% stocashtic parameters

rhoint=0.8;
rhogov=0.8;
rhozy=0.9;
rhozi=0.5;
rhozmk=0.5;
rhozwmk=0.5;
rhozrp=0.9;


