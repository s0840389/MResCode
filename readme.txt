

Folder structure

    1) charts 

    Contains matlab and R codes to generate figures in text

        charts.m - creates most IRF charts in section 4

        decomp_charts.m - Does the partial eq consumption decomp of the HANK IRF's (figure 7)

        charts.r - creates consumption decomp chart (figure 7)

        VAR/main.m - Estimates the VAR in figure 1. [Data included] 

        ls_y_correlation.r - correlation between output and labour share 

                

    2) HANK 

    Contains 3 HANK models solved using the peturbation algorithms Bayer, Born & Luetticke

    NK - medium scale new keynesian model
    NKYN - model augmented with expansionary labour but flexible labour markets Wy=We
    NKYN2W - model augmented with expansionary labour but inflexible labour markets Wy~=We

steady states are saved in directory '../../../steadysteates' new users need to resolve the steady state and save in a chosen directory. 


    3) RANK

     main.m solves 4 models and produces IRF's for those models using a first order peturbation

        NK -  a text book medium scale model

        NKYN - model augmented with expansionary labour

        NKcap - two agent baby hank model with workers and capatalists

        NKYNcap - two agent baby hank model with workers and capatalists augmented with expansionary labour


    4) DynareCode 
    
    Contains the dynare mod files and data used to estimate the RANK parameters. 

    DyamicsNK.mod - main file for the NK model
    DyamicsNKYN.mod - main file for the NK-YN model
