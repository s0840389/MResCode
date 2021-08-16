%% Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Household Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.beta        = 0.98;     % Discount factor
par.xi          = 4;          % CRRA
par.gamma       = 2.00;          % Inverse Frisch elasticity
par.nu          = 0.065;          % Prob. of trade given adj. decision

% Income Process
par.rhoH        = 0.98;    % Persistence of productivity
par.sigmaH      = 0.06;    % STD of productivity shocks
mpar.in         = 0.00063;  % Prob. to become superstar/entrepreneur
mpar.out        = 0.0625;   % Prob. to become worker again


par.SSsh=0.12;             % super star share of labour income

par.lumpshare=0.0; % share of dividends that are lump sum versus goes to investment account

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm Side Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.PIdshr=0.05;
par.Labshr=0.62;
par.shareE=0.2;


par.eta         = 1/(par.PIdshr+par.Labshr*par.shareE);

par.mu          = (par.eta-1)/par.eta;       % Markup
%par.alpha       = 2/3/par.mu;  % Labor share
par.delta       = 0.07/4;     % Depreciation rate
par.phi         = 0.23;        % Capital adj costs

par.thetay     = 1.00;     % scale Y



 lyle_ratio=par.shareE/(1-par.shareE);
 
 scaley=1/par.mu*par.Labshr/(1+lyle_ratio); % thetay*(1-alphay)
 
 scalen=1/(1-par.mu)*par.Labshr/(1+1/lyle_ratio); % thetan*(1-alphan)
    
 par.thetay=1;
 
 par.alphay=1-scaley;
 
 par.thetan=scalen;

% Tax Schedule

par.GY=0.18;

par.tau         = 1-(par.GY+1.6*0.025/4)/par.Labshr;   % Proportional tax on labor and profit income 

par.rho_tax=0.55;
par.gamma_taxB=0.777;
par.gamma_taxY=2.646;

% Central Bank Policy
par.theta_pi    = 1.8;  % Reaction to inflation
par.rho_R       = 0.8;  % Inertia
par.theta_y=0.1;        %reaction to output

% Debt rule [turned off]
par.gamma_pi    = -1.13;   % Reaction to inflation
par.gamma_Y     =-0.716; % Reaction to tax revenue
par.gamma_B       = 0.14;  % Autocorrelation


% Phillips Curve
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)

% wage philips curve

par.epsilon_w=11;

par.Nw=4; % average wage duration

par.thetaw=(par.Nw-1)/(par.Nw); % calvo wage param

par.phiw=par.thetaw*(par.epsilon_w-1)*par.tau/((1-par.thetaw)*(1-par.beta*par.thetaw)); % wage adjustmetn cost

par.kappa_w=(1-par.thetaw)*(1-par.beta*par.thetaw)/par.thetaw;


%% Returns
par.PI  = 1.00^.25;     % Gross inflation
par.RB  = (par.PI*1.025)^0.25;    % Real return times inflation

par.ABS = 0;                    % Loan to value ratio max.
par.borrwedge = par.PI*(1.1^0.25-1); % Wedge on borrowing beyond secured borrowing

par.Q  = 1;

%% Grids
% Idiosyncratic States
mpar.nm         = 60;
mpar.nk         = 60;
mpar.nh         = 15;
mpar.tauchen    ='importance';


%% Numerical Parameters
mpar.crit    = 1e-10;

