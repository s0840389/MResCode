clear all
close all

NKcol =  [0.2 0.2 1.0];
YNcol =  [1 0 0.83];
YN2Wcol= [0.0 0.8 0.0];

%% IRFs


IRFpos=[4,3,2]';
IRFname={'G','TFP','IST'};
IRFname2={'Government spending','TFP','Investment-specific technology'};


% NK model
load('../../../steadystates/NKfund_60_15.mat','hx','gx','mpar')

mparNK=mpar;
hxNK=hx;
gxNK=gx;

% YN model
load('../../../steadystates/YNfund_60_15.mat','hx','gx','par','mpar')

for i=1:size(IRFpos,1)

x0=zeros(mpar.numstates,1);
x0(end-IRFpos(i))=-0.01;

MX=[eye(length(x0));gxNK];
IRFXX=[];
x=x0;

% y=ybar+gx*(x-xbar)
% x'=xbar+hx(x-xbar)

for t=1:mpar.maxlag
    IRFXX(:,t)=(MX*x)';
    x=hxNK*x*exp(-max(0,t-25)*1/250);
end

eval(sprintf('IRF_NK_%s=IRFXX;',string(IRFname(i))))


end




for i=1:size(IRFpos,1)

x0=zeros(mpar.numstates,1);
x0(end-IRFpos(i))=-0.01;

MX=[eye(length(x0));gx];
IRFXX=[];
x=x0;

% y=ybar+gx*(x-xbar)
% x'=xbar+hx(x-xbar)

for t=1:mpar.maxlag
    IRFXX(:,t)=(MX*x)';
    x=hx*x*exp(-max(0,t-25)*1/250);
end

eval(sprintf('IRF_YN_%s=IRFXX;',string(IRFname(i))))


end

%% figure 

figure(1)
clf

% Output

  for i=1:size(IRFpos,1)
      
    subplot(4,3,i)

    hold on
    
    eval(sprintf('plot(100*IRF_NK_%s(end-mparNK.oc+3,1:end-1),"Color",NKcol,"LineWidth",1.8)',string(IRFname(i))))

    eval(sprintf('plot(100*IRF_YN_%s(end-mpar.oc+3,1:end-1),"Color",YNcol,"LineWidth",1.8)',string(IRFname(i))))

    hline=refline(0,0);
    hline.Color='black';

    
    if i==1
    ylabel('Output')
    end
    
    title(string(IRFname2(i)))

  end

% Inflation

  for i=1:size(IRFpos,1)
      
    subplot(4,3,3+i)

    hold on
    
    eval(sprintf('plot(100*IRF_NK_%s(end-mparNK.oc+2,1:end-1),"Color",NKcol,"LineWidth",1.8)',string(IRFname(i))))

    eval(sprintf('plot(100*IRF_YN_%s(end-mpar.oc+2,1:end-1),"Color",YNcol,"LineWidth",1.8)',string(IRFname(i))))

    hline=refline(0,0);
    hline.Color='black';

    
    if i==1
    ylabel('Inflation')
    end
    

  end


% Labour share

  for i=1:size(IRFpos,1)
      
    subplot(4,3,6+i)

    hold on
    
    eval(sprintf('LS=100*par.Labshr*(IRF_NK_%s(end-mparNK.oc+7,1:end-1)+IRF_NK_%s(mparNK.numstates-5,2:end)-IRF_NK_%s(end-mparNK.oc+3,1:end-1));',string(IRFname(i)),string(IRFname(i)),string(IRFname(i))))
    
    plot(LS,"Color",NKcol,"LineWidth",1.8)
    
 eval(sprintf('LS=100*par.Labshr*(IRF_YN_%s(end-mpar.oc+7,1:end-1)+IRF_YN_%s(mpar.numstates-5,2:end)-IRF_YN_%s(end-mpar.oc+3,1:end-1));',string(IRFname(i)),string(IRFname(i)),string(IRFname(i))))
    
    plot(LS,"Color",YNcol,"LineWidth",1.8)
    
    hline=refline(0,0);
    hline.Color='black';

    
    if i==1
    ylabel('Labour share')
    end
    

  end

  
% Productivity

  for i=1:size(IRFpos,1)
      
    subplot(4,3,9+i)

    hold on
    
    eval(sprintf('LS=100*(-IRF_NK_%s(end-mparNK.oc+7,1:end-1)+IRF_NK_%s(end-mparNK.oc+3,1:end-1));',string(IRFname(i)),string(IRFname(i))))
    
    plot(LS,"Color",NKcol,"LineWidth",1.8)
    
 eval(sprintf('LS=100*(-IRF_YN_%s(end-mpar.oc+7,1:end-1)+IRF_YN_%s(end-mpar.oc+3,1:end-1));',string(IRFname(i)),string(IRFname(i))))
    
    plot(LS,"Color",YNcol,"LineWidth",1.8)
    
    hline=refline(0,0);
    hline.Color='black';

    
    if i==1
    ylabel('Productivity')
    legend('HANK','HANK-YN','Location','southeast')
    end
    

  end

  
%%  figure(1)

h = gcf;
set(h,'Units','Inches');
h.Position(1) =h.Position(1)*0;
h.Position(2) =h.Position(2)*0;
h.Position(3) =h.Position(3)*1.4;
h.Position(4) =h.Position(4)*1.4;
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'IRF_othershocks','-dpdf','-r0')