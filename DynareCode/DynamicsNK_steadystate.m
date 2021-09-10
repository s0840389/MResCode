
function [ys,check] = DynamicsNK_steadystate(ys,exe);
global M_

%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;
%%
%% END OF THE FIRST MODEL INDEPENDENT BLOCK.


pkl=1/bbeta+ddelta-1-0.01;
pkh=0.15;

diff=10;
tol1=10e-8;
iter=0;

while abs(diff)>tol1

iter=iter+1;

% (1) pk 

	pk=0.5*(pkl+pkh);

% (2) inflation

   pit=0;
   
% (3) interest rate

   int=0; %(16)

% (4) Production labour

	ly=log(1);

% (5) marginal cost

	mc=log(1/(mu_p-(pit-bbeta*pit)/kappa_p));
    
% (6) technology

	zy=0;

% (7) capital

	k=log((pk/(exp(mc)*tthetay*aalpha*exp(zy)*exp(ly)^((1-aalpha)*tthetay)))^(1/(aalpha*tthetay-1)));

% (8) firm output

	y=log(((exp(k))^aalpha*exp(ly)^(1-aalpha))^tthetay); % firm level production

% (9) wages
 
	w=log(exp(mc)*tthetay*(1-aalpha)*exp(y)/exp(ly)); 
 
  % (12) dividend
  Inv=log((ddelta+ctrend)*exp(k));
  
 PId=log(exp(y)-exp(w)*exp(ly)-pk*exp(k));
 
% (13) no arbritrage


qk=0;

	ra=(pk+exp(qk)*(1-ddelta)) / exp(qk) -1;

% (14) stock price

	q=log(exp(PId)/ra);

% (15) fund value

	a=log(exp(q)+exp(k)*exp(qk));

% (16) deposit

	d=(ctrend-ra)*exp(a);

	B=log(BY*exp(y));

	taxrate=taxss;

	taxrev=exp(w)*exp(ly)*(1-taxrate);

	g=log(taxrev+exp(B)*(1+ctrend -(1+int)/(1+pit)));

% (17) consumption

	adjcost=cchi0*abs(d) + cchi1*abs(d)^cchi2*exp(a)^(1-cchi2);

c=log(exp(y)-exp(g)-exp(Inv));

%	c=log(exp(w)*exp(ly)*taxrate+(1+ctrend)*exp(B)*int/(1+pit)-d);

 % other variables
 
	sy=exp(ly)*exp(w)/(exp(y));


yf=y;

zi=0;
zrp=0;
zmk=0;
zwmk=0;


pitw=0;

dy=ctrend*100;
dc=ctrend*100;
dinve=ctrend*100;
dw=ctrend*100;
pinfobs = pitstar*100;
robs =    intstar*100;
labobs = constelab; 

% difference

	diff= (1+ra) - (1+ctrend)^(ssigma)/bbeta - euler2(d,cchi1,cchi2,cchi0,a);

if diff>0

	pkh=pk;

else

	pkl=pk;

end

if iter>200
break
end

end
 
lss=ly;
wss=w;
gss=g;
Bss=B;
yss=y;
mcss=mc;


yf=y;
cf=c;
lyf=ly;
wf=w;
pkf=pk;
qkf=qk;
kf=k;
Invf=Inv;
%%
%% END OF THE MODEL SPECIFIC BLOCK.


%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state values of the endogenous variables of
%% the model.
%%

NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
for i = 1:NumberOfEndogenousVariables                         % Loop...
  varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
  eval(['ys(' int2str(i) ') = ' varname ';']);               %    Get the steady state value of this variable.
 % varname 
 % [ys(i) exp(ys(i))]

end                                                           % End of the loop.

% update parameters

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end


end

function x = euler2(d,cchi1,cchi2,cchi0,a)
	
	xd=cchi0+cchi2*cchi1*abs(d)^(cchi2-1)*exp(a)^(1-cchi2);
    xa=(1-cchi2)*cchi1*abs(d)^cchi2*exp(a)^(-cchi2);
	x=xa/(1+xd);
    end
%%
%% END OF THE SECOND MODEL INDEPENDENT BLOCK.