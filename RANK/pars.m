

%% households

p.beta=0.9938; % discount
p.sigma=4.0;	% crra
p.psi=2.0;	% labour elast (inverse frish)


p.epsilon_w=11;	% epsilon wagws
p.mu_w=p.epsilon_w/(p.epsilon_w-1); % wage markup

p.shareE=0.2; % share of expansionary workers

p.shrCap=0.1;
p.betaCap=0.98; % discount

%% firms

p.labshr=0.62;
p.profshr=0.05;


p.epsilon_p=1/p.profshr;	% epsilon prices
p.mu_p=p.epsilon_p/(p.epsilon_p-1); % price markup

p.alpha=1-p.labshr*p.mu_p;	% capital elasticity
p.delta=0.07/4;
p.thetay=1.00;	% scale Y
p.thetan=0.90; % scale N

%% government

p.phipi=1.8; % taylor rule

p.phiy=0.1;

p.rhoint=0.8; % smoothing

p.pitstar=0.00/4; % inflation target

p.intstar=(1/p.beta)*(1+p.pitstar)-1; % r*


% Tax Schedule


p.BY=1.6; %debt to gdp

p.taul=1-(0.175+p.BY*p.intstar)/p.labshr;


% tax rules
p.rho_tax=0.55; %autocorrelation
p.gamma_taxB=0.777; %debt
p.gamma_taxY=2.646; %output

% Debt rule
p.gamma_pi    = -1.13;   % Reaction to inflation
p.gamma_Y     =-0.716; % Reaction to tax revenue
p.gamma_B       = 0.14;  % Autocorrelation


%% frictions

p.tau=1.8; % investment adjustment costs [3 percent on impact] 

p.Np=4; % average price duration

p.theta=(p.Np-1)/p.Np; % calvo param


p.phi=p.theta*(p.epsilon_p-1)/((1-p.theta)*(1-p.beta*p.theta)); % price adjustment cost

p.kappa_p=(1-p.theta)*(1-p.beta*p.theta)/p.theta;

p.Nw=4; % average wage duration

p.thetaw=(p.Nw-1)/p.Nw; % calvo wage param

p.phiw=p.thetaw*(p.epsilon_w-1)*p.taul/((1-p.thetaw)*(1-p.beta*p.thetaw)); % wage adjustmetn cost

p.kappa_w=(1-p.thetaw)*(1-p.beta*p.thetaw)/p.thetaw;

p.chi0=0*0.05;
p.chi1=0*0.8;
p.chi2=1.4;

%% stocashtic parameters

p.rhozy=0.9;

p.se_int=0.01;
p.se_zy=0.01;

%% other parameters

p.Nss=1; % steady state labour
p.Eshr=0.2; % share of expansionary labour
