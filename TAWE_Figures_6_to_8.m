%% Plots for TAWE tremor suppression with PD, MPC (with and without BMFLC)(Figs 6,7,8)

% Simulation code for "Model Predictive Control for Tremor Suppressing Exoskeleton"
% by Subham Samal (subhamsamal@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech

% The authors thank Jiamin Wang (jmechw@vt.edu) for their constant guidance
% and support throughout the project

% This script contains the code for plotting Figs 6, 7 and 8 for the
% wearable wrist exoskeleton simulations for tremor suppression 
PathSetup;
%% Load Data and Parameters

drawingParam;

% load data for PD Controller Simulations
zero_PD = load('zero_track_TAWE_PD.mat');
zero_pd_tr = zero_PD.simStateData;

% load data for MPC Controller without BMFLC
zero_mpc_tr = load('zero_track_TAWE_MPC.mat');
zero_mpc_inp = load('zero_track_TAWE_MPC_inp.mat');

zero_mpc_tr = zero_mpc_tr.simStateData;
zero_mpc_inp = zero_mpc_inp.simInputData;

% load data for MPC with BMFLC
zero_MPC_BMFLC = load('zero_track_bmflc.mat');

zero_mpc_bmflc_tr = zero_MPC_BMFLC.simStateData;
hh = 0.002;
tSpan = 0:hh:60;

%% PD Error (Fig 6)

% Generate plots
Fig=figure(6);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,1,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zeros(1,30000)-zero_pd_tr(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 10])
ylim([-0.15 0.15])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$\epsilon_{\theta,\mathrm{RUD}}$$','$$\epsilon_{\theta,\mathrm{RUD}}$$$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(2,1,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zeros(1,30000)-zero_pd_tr(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 10])
ylim([-0.15 0.15])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$\epsilon_{\theta,\mathrm{FE}}$$','$$\epsilon_{\theta,\mathrm{FE}}$$$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;


%% PD Error vs MPC Fig 7

Fig=figure(7);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zeros(1,30000)-zero_pd_tr(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(4),'Color',lineColorDefault{2})
hold on
plot(tSpan(1:30000),zeros(1,30000)-zero_mpc_tr(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(4),'Color',lineColorDefault{3})
xlim([0 10])
ylim([-0.15 0.15])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon_{RUD}}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('PD','MPC','interpreter','latex','FontSize',15,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zeros(1,30000)-zero_pd_tr(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(4),'Color',lineColorDefault{2})
hold on
plot(tSpan(1:30000),zeros(1,30000)-zero_mpc_tr(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(4),'Color',lineColorDefault{3})
xlim([0 10])
ylim([-0.12 0.2])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon_{FE}}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('PD','MPC','interpreter','latex','FontSize',15,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=3;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zero_mpc_inp(3,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(4),'Color',lineColorDefault{4})
xlim([0 10])
ylim([-2 2])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$u_\mathrm{RUD}$$ (Nm)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('MPC','interpreter','latex','FontSize',15,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=4;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zero_mpc_inp(4,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(4),'Color',lineColorDefault{4})
xlim([0 10])
ylim([-2 2])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$u_\mathrm{FE}$$ (Nm)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('MPC','interpreter','latex','FontSize',15,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

%% With and without BMFLC Fig 8

Fig=figure(8);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,1,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zeros(1,30000)-zero_mpc_tr(1,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan(1:30000),zeros(1,30000)-zero_mpc_bmflc_tr(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 10])
ylim([-0.1 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon_{RUD}}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('Without BMFLC','With BMFLC','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(2,1,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),zeros(1,30000)-zero_mpc_tr(2,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan(1:30000),zeros(1,30000)-zero_mpc_bmflc_tr(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 10])
ylim([-0.1 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon_{FE}}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('Without BMFLC','With BMFLC','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;


