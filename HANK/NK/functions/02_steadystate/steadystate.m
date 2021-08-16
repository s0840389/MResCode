function [c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr,R_fc,W_fc,Profits_fc,Output,grid,par,qs]=steadystate(P_H,grid,meshes,mpar,par)
%     Prepare items for EGM
%     Layout of matrices:
%       Dimension 1: money m
%       Dimension 2: capital k
%       Dimension 3: stochastic human capital h

% 1) Construct relevant return matrices
joint_distr = ones(mpar.nm,mpar.nk,mpar.nh)/(mpar.nh*mpar.nk*mpar.nm);

[ ~,~,~,Profits_fc,WW,RR,RBRB,~] = factor_returns(meshes,grid,par,mpar);

% 2) Guess initial policies
[c_a_guess,c_n_guess, psi_guess,~]=policyguess(meshes,WW,RR,RBRB,par,mpar,Profits_fc,grid.K);

KL=0.7*grid.K;
KH=1.3*grid.K;

disp('Calculating excessL')
[excessL]=excessK(KL,c_a_guess,c_n_guess,psi_guess,joint_distr, grid,P_H,mpar,par,meshes);
excessL
disp('calculating excessH')
[excessH]=excessK(KH,c_a_guess,c_n_guess,psi_guess,joint_distr, grid,P_H,mpar,par,meshes);
excessH

if sign(excessL)==sign(excessH)
    disp('ERROR! Sign not diff')
end

%% Brent
fa=excessL; fb=excessH;
a=KL;
b=KH;

if fa*fb>0
    error('f(a) und f(b) sollten unterschiedliche Vorzeichen haben');
end

c=a; fc=fa;   %At the beginning: c = a

c=a; fc=fa; d=b-a; e=d;

iter=0;
maxiter=1000;

while iter<maxiter
    iter=iter+1
    
    
    if fb*fc>0
        c=a; fc=fa; d=b-a; e=d;
    end
    
    if abs(fc)<abs(fb)
        a=b; b=c; c=a;
        fa=fb; fb=fc; fc=fa;
    end
    
    tol=2*eps*abs(b)+mpar.crit; m=(c-b)/2; % tolerance
    
    if (abs(m)>tol) && (abs(fb)>0)
        
        if (abs(e)<tol) || (abs(fa)<=abs(fb))
            d=m; e=m;
        else
            s=fb/fa;
            if a==c
                p=2*m*s; q=1-s;
            else
                q=fa/fc; r=fb/fc;
                p=s*(2*m*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            end
            if p>0
                q=-q;
            else
                p=-p;
            end
            s=e; e=d;
            if ( 2*p<3*m*q-abs(tol*q) ) && (p<abs(s*q/2))
                d=p/q;
            else
                d=m; e=m;
            end
        end
        
        a=b; fa=fb;
        
        if abs(d)>tol
            b=b+d;
        else
            if m>0
                b=b+tol;
            else
                b=b-tol;
            end
        end
    else
        break;
    end
    
    [fb,c_n_guess,~,c_a_guess,~,~,psi_guess,joint_distr,~]=excessK(b,c_a_guess,c_n_guess,psi_guess,joint_distr, grid,P_H,mpar,par,meshes);
    fb
end

Kcand=b;

%% Update

[excess,c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr, R_fc,W_fc,Profits_fc,Output,N,qs,par]=excessK(Kcand,c_a_guess,c_n_guess,psi_guess,joint_distr, grid,P_H,mpar,par,meshes);

grid.N=N;
grid.K=Kcand;
end

function [excess,c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,joint_distr, R_fc,W_fc,Profits_fc,Output,N,qs,par,ly,le,Mg]=excessK(K,c_a_guess,c_n_guess,psi_guess,joint_distr, grid,P_H,mpar,par,meshes)

grid.K=K;

[ N,R_fc,W_fc,Profits_fc,WW,RR,RBRB,Output,qs,par] = factor_returns(meshes,grid,par,mpar);



[~,~, ~,inc]=policyguess(meshes,WW,RR,RBRB,par,mpar,Profits_fc,grid.K);

% Solve Policies and Joint Distribution
disp('Solving household problem by EGM')
tic
[c_n_guess,m_n_star,c_a_guess,m_a_star,cap_a_star,psi_guess,distPOL]=...
    policies_SS(c_a_guess,c_n_guess,psi_guess, grid, inc, RR,RBRB,P_H,mpar,par,meshes);
disp(distPOL)
toc

disp('Calc Joint Distr')
tic
[joint_distr,distJD]=JDiteration(joint_distr,m_n_star,m_a_star,cap_a_star,P_H,par,mpar,grid);
joint_distr=reshape(joint_distr,[mpar.nm mpar.nk mpar.nh]);
disp(distJD)
toc
AggregateCapitalDemand=sum(grid.k.*sum(sum(joint_distr,1),3));
excess=(grid.K+qs-AggregateCapitalDemand);
end
