close all

%% Plotting

figurename=['IRF_Y_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_Y,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);
printpdf(gcf,['latex/' figurename])

figurename=['IRF_C_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_C,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_I_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_I,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_G_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_G,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_Deficit_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,Deficit,'b-','LineWidth',4.5)
ylabel('Percentage Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_K_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_K,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_M_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_M,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_H_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_H,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_S_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_S,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])
%%
figurename=['IRF_RBPI_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_RB,'--','LineWidth',4.5)
hold on
plot(1:mpar.maxlag-1,IRF_RBREAL,'b-','LineWidth',4.5)
hold on
legend({'nominal','real'},'Location','NorthEast','AutoUpdate','off')
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_RB_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_RB,'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_PI_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_PI,'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_Q_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_Q,'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_D_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_D,'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

figurename=['IRF_LP_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-2,IRF_LP,'b-','LineWidth',4.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])


figurename=['IRF_N_' casename '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1000 1000])
plot(1:mpar.maxlag-1,IRF_N,'b-','LineWidth',4.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',45)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',45)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',45);  
printpdf(gcf,['latex/' figurename])

