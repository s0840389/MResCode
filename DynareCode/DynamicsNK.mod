// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating conditional expectation with 
// polynomials)
//


//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parametersNK.mod"

//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------



@#include "variablesNK.mod"

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equationsNK.mod"

end;

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;

var eint;  stderr 0.01;
var egov;  stderr 0.01;
var ezy;  stderr 0.01;
var ezi;  stderr 0.01;
var ezmk;  stderr 0.005;
var ezwmk;  stderr 0.005;
var ezrp;  stderr 0.005;
end;

options_.steadystate.nocheck = 1;
steady;

// Compute steady state (nocheck option ensures that Dynare runs even if steady
// state only computed approximately, i.e., with small numerical error)
//steady(nocheck);

// Check regularity conditions (turn on to check)
//check;
//model_diagnostics;
//model_info;

estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF

  
	// steady state params 

	ctrend,0.0042,0.001,0.008,NORMAL_PDF,0.0038,0.005;
	pitstar,0.0075,0.001,0.02,GAMMA_PDF,0.00625,0.001;
	constelab,1.2918,-10.0,10.0,NORMAL_PDF,0.0,2.0;

	cbeta,0.0074,0.0001,0.02,GAMMA_PDF,0.0025,0.001;
	ssigma, 2,0.25,4,NORMAL_PDF,1.5,1;
	ppsi,2,0.1,6,NORMAL_PDF,2.00,1;
	
	// dynamics params

	kappa_p,0.10,0.01,0.5,GAMMA_PDF,0.10,0.02;
	kappa_w,0.10,0.01,0.5,GAMMA_PDF,0.10,0.02;
	ttau,2.0,0.01,10,GAMMA_PDF,2,2;

	phipi,1.5,1.05,4,NORMAL_PDF,1.5,0.3;
	phiy,0.1,0.01,1,NORMAL_PDF,0.1,0.05;
	phiyy,0.0893,0.001,0.5,NORMAL_PDF,0.125,0.05;

	//rho_tax,.5264 ,.01,.9999,BETA_PDF,0.5,0.20;
	//gamma_taxB,0.1710,0.1,2,NORMAL_PDF,0.5,1;
	//gamma_taxY,0.15,0.1,3,NORMAL_PDF,1.0,2;
	
	// shocks 

	stderr eint,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhoint,.5 ,.01,.90,BETA_PDF,0.5,0.15;

	stderr egov,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhogov,.5 ,.01,.9999,BETA_PDF,0.5,0.20;

	stderr ezy,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhozy,.5 ,.01,.9999,BETA_PDF,0.5,0.20;

	stderr ezi,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhozi,.5 ,.01,.9999,BETA_PDF,0.5,0.20;

	stderr ezmk,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhozmk,.5 ,.01,.9999,BETA_PDF,0.5,0.20;

	stderr ezwmk,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhozwmk,.50 ,.01,.9999,BETA_PDF,0.5,0.20;

	stderr ezrp,0.01,0.0001,0.03,INV_GAMMA_PDF,0.01,0.02;
	rhozrp,.5 ,.01,.9999,BETA_PDF,0.5,0.20;


end;

varobs dy  pinfobs robs dc dinve dw labobs;


// Find mode
//estimation(datafile=usmodel_data,mode_compute=1,mode_check,first_obs=72,presample=76,prefilter=0,mh_replic=10000,mh_nblocks=2,mh_jscale=0.3,mh_drop=0.2);

// MH
//estimation(datafile=usmodel_data,mode_compute=0,mode_File=DynamicsNK_mode,first_obs=72,presample=76,prefilter=0,mh_replic=75000,mh_nblocks=2,mh_jscale=0.3,mh_drop=0.2);

// IRFs
//estimation(datafile=usmodel_data,mode_compute=0,mode_File=DynamicsNK_mode,load_mh_file,first_obs=72,presample=76,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.3,mh_drop=0.2,nodiagnostic,plot_priors=0,bayesian_irf);

//Var decomp
//estimation(datafile=usmodel_data,mode_compute=0,mode_File=DynamicsNK_mode,load_mh_file,first_obs=72,presample=76,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.3,mh_drop=0.2,nodiagnostic,plot_priors=0,moments_varendo,conditional_variance_decomposition=[1 4 8 12]);

// simulation
//estimation(datafile=usmodel_data,mode_compute=0,mode_File=DynamicsNK_mode,load_mh_file,first_obs=72,presample=76,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.3,mh_drop=0.2,nodiagnostic,plot_priors=0);

//stoch_simul(periods=1000,order=1,hp_filter=1600) y sy;


estimation(datafile=usmodel_data,mode_compute=0,mode_File=DynamicsNK_mode,load_mh_file,first_obs=72,presample=76,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.3,mh_drop=0.2,plot_priors=0);
