

load('../../../steadystates/NKfund_60_15.mat')


usgdp=20*10e12/10e6;
uswf=157;

scal=usgdp/uswf;

figure(1)

subplot(2,1,1)

   h = area(grid.k/Output,capdist,'FaceAlpha',0.20);
            
   xlim([-5 inf])

title('Illiquid wealth distribution')
   
ylabel('PDF')
xlabel('Wealth relative to Output per capita')

subplot(2,1,2)
h = area(grid.m/Output,moneydist,'FaceAlpha',0.20);
            
title('Liquid wealth distribution')

ylabel('PDF')
xlabel('Wealth relative to Output per capita')
 xlim([-5 inf])

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
%h.Position(3) =h.Position(3)*1.2;
%h.Position(4) =h.Position(4)*0.9;
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'A2_ssdist','-dpdf','-r0')

%%

Income=meshes.k*par.R+meshes.m*(par.RB-1)+par.W*meshes.h*par.tau;

[income_sort, IX] = sort(Income(:));
income_pdf        = joint_distr(IX);


figure(2)

   h = area(income_sort/Output,income_pdf,'FaceAlpha',0.20);
            
  
