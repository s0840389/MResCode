function [P_H,grid,par]=stochastics_variance(par, mpar, grid)
% stochastics_variance generates transition probabilities for:
% h: P_H

% First for human capital (prodcution)
% Generate transition probabilities and grid
[hgridy,P_Hy,boundsHy] = Tauchen(par.rhoHy,mpar.nhy-1,1, 0, mpar.tauchen); % LR variance = 1
% Correct long run variance for *human capital*
hgridy               = hgridy*par.sigmaHy/sqrt(1-par.rhoHy^2);
hgridy               = exp(hgridy); % Levels instead of Logs

P_Hy     = transition(mpar.nhe-1,par.rhoHy,sqrt(1-par.rhoHe^2),boundsHy);

% solve for superstar income

P_Hy_erg=[1 zeros(1,size(hgridy,2)-1)]*P_Hy^1000;

hSS=par.SSsh*sum(P_Hy_erg.*hgridy)/(0.01-par.SSsh*0.01);

% Transitions to superstar/entrepreneur state
P_Hy=[P_Hy repmat(mpar.in,[mpar.nhy-1 1])];
lastrow=[repmat(0,[1, mpar.nhy-1]) 1-mpar.out];
lastrow(ceil(mpar.nhy/2))=mpar.out;
P_Hy=[P_Hy; lastrow];
P_Hy=P_Hy./repmat(sum(P_Hy,2),[1 mpar.nhe]);

hgridy=[hgridy hSS];

% First for human capital (expansion)
% Generate transition probabilities and grid
[hgride,P_He,boundsHe] = Tauchen(par.rhoHe,mpar.nhe-1,1, 0, mpar.tauchen); % LR variance = 1
% Correct long run variance for *human capital*
hgride               = hgride*par.sigmaHe/sqrt(1-par.rhoHe^2);
hgride               = exp(hgride); % Levels instead of Logs

P_He     = transition(mpar.nhe-1,par.rhoHe,sqrt(1-par.rhoHe^2),boundsHe);

% Transitions to superstar/entrepreneur state
P_He=[P_He repmat(mpar.in,[mpar.nhe-1 1])];
lastrow=[repmat(0,[1, mpar.nhe-1]) 1-mpar.out];
lastrow(ceil(mpar.nhe/2))=mpar.out;
P_He=[P_He; lastrow];
P_He=P_He./repmat(sum(P_He,2),[1 mpar.nhe]);

hgride=[hgride hSS];

grid.h=[hgridy hgride];

P_H=[P_Hy zeros(mpar.nhy,mpar.nhe); 
    zeros(mpar.nhe,mpar.nhy) P_He];

Pauxy=P_Hy^1000^1000;
hh=Pauxy*hgridy';
par.Hy=hh(1); % Total prod Employment

Pauxe=P_He^1000^1000;
hh=Pauxe*hgride';
par.He=hh(1); % Total expansion Employment

%par.profitshare=Paux(end,end)^-1; % Profit per household
grid.boundsH=boundsHy;

end