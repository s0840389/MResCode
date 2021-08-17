

Folder structure

    1) charts 

    Contains matlab and R codes to generate figures in text

    charts.m - creates charts in section 4

    decomp_charts.m - Does the partial eq consumption decomp of the HANK IRF's

    charts.r - creates consumption decomp chart

    figure1.m - constructs VAR IRF for the labour share [ need to add code for actual VAR]

    2) HANK 

    Contains 3 HANK models solved using the peturbation algorithms Bayer, Born & Luetticke

    NK - medium scale new keynesian model
    NKYN - model augmented with expansionary labour but flexible labour markets Wy=We
    NKYN2W - model augmented with expansionary labour but inflexible labour markets Wy~=We

steady states are saved in directory '../../../steadysteates' new users need to resolve the steady state and save in a chosen directory. 


    3) RANK

     main.m solves 4 models and produces IRF's for those models

    NK -  a text book medium scale model

    NKYN - model augmented with expansionar labour

    NKcap - two agent baby hank model with workers and capatalists

    NKYNcap - two agent baby hank model with workers and capatalists augmented with expansionary labour

