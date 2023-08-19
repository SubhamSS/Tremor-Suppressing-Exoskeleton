%% Plots for Arm and TAWE trajectory tracking with MPC, no tremor (Figs 3,4,5)

% Simulation code for "Model Predictive Control for Tremor Suppressing Exoskeleton"
% by Subham Samal (subhamsamal@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech

% The authors thank Jiamin Wang (jmechw@vt.edu) for their constant guidance
% and support throughout the project

% This script contains the code for plotting Figs 3, 4 and 5 for the
% wearable wrist exoskeleton and arm simulations for trajectory tracking using MPC
PathSetup;
%% Load Data and Parameters
drawingParam;

% load reference data
load('simTAWERef','refData','refdtData','refddtData');

% load data for arm simulation (Fig 3)
arm_sim_data = load('arm_no_trem.mat');
arm_sim_states = arm_sim_data.simStateData;

% load data for TAWE Simulation
tawe_sim_data = load('tawe_traj_track_no_trem.mat');
tawe_inp_data = load('tawe_traj_track_no_trem_inp.mat');

tawe_sim_states = tawe_sim_data.simStateData;
tawe_sim_inp = tawe_inp_data.simInputData;
hh = 0.002;
tSpan = 0:hh:60;

%% Figure 3 Arm trajectory tracking

Fig=figure(3);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),refData(1,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan(1:30000),arm_sim_states(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-1.1 1.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\theta_\mathrm{RUD}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$r_{\theta,\mathrm{RUD}}$$','$$\theta_\mathrm{RUD}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),refData(2,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan(1:30000),arm_sim_states(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-1.1 1.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\theta_\mathrm{FE}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$r_{\theta,\mathrm{FE}}$$','$$\theta_\mathrm{FE}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=3;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),refData(1,1:30000)-arm_sim_states(1,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
xlim([0 60])
ylim([-0.1 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylh = ylabel(strcat('\textbf{$$\mathbf{\epsilon}$$  (rad)}'),'interpreter','latex');
ylh.Position(1) = ylh.Position(1) + 0.1;
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$\epsilon_{\theta,\mathrm{RUD}}$$','$$\epsilon_{\theta,\mathrm{FE}}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=4;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),refData(2,1:30000)-arm_sim_states(2,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
xlim([0 60])
ylim([-0.1 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylh = ylabel(strcat('\textbf{$$\mathbf{\epsilon}$$  (rad)}'),'interpreter','latex');
ylh.Position(1) = ylh.Position(1) + 0.1;
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$\epsilon_{\theta,\mathrm{FE}}$$','$$\epsilon_{\theta,\mathrm{FE}}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

%% Figure 4 Arm and TAWE-ForeArm Comparison plot

Fig=figure(4);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),arm_sim_states(1,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{4})
hold on
plot(tSpan(1:30000),tawe_sim_states(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-1.1 1.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\theta_\mathrm{RUD}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',20,'FontWeight','bold')
legend('Arm','TAWE-Arm','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),arm_sim_states(2,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{4})
hold on
plot(tSpan(1:30000),tawe_sim_states(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-1.6 1.4])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\theta_\mathrm{FE}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',20,'FontWeight','bold')
legend('Arm','TAWE-Arm','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=3;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),refData(1,1:30000)-arm_sim_states(1,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{4})
hold on
plot(tSpan(1:30000),refData(1,1:30000)-tawe_sim_states(1,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-0.05 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylh = ylabel(strcat('\textbf{$$\mathbf{\epsilon_{RUD}}$$  (rad)}'),'interpreter','latex');
ylh.Position(1) = ylh.Position(1) - 0.05;
set(gca,'fontname','times new roman','FontSize',20,'FontWeight','bold')
legend('Arm','TAWE-Arm','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=4;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),refData(2,1:30000)-arm_sim_states(2,1:30000),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{4})
hold on
plot(tSpan(1:30000),refData(2,1:30000)-tawe_sim_states(2,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-0.05 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylh = ylabel(strcat('\textbf{$$\mathbf{\epsilon_{FE}}$$  (rad)}'),'interpreter','latex');
ylh.Position(1) = ylh.Position(1) - 0.05;
set(gca,'fontname','times new roman','FontSize',20,'FontWeight','bold')
legend('Arm','TAWE-Arm','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

%% Figure 5 Tawe Inputs

Fig=figure(5);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,1,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),tawe_sim_inp(3,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{4})
xlim([0 60])
ylim([-3 3])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$u_\mathrm{RUD}$$ (Nm)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(2,1,jj,0.3,0.15,0.05,0.05));
plot(tSpan(1:30000),tawe_sim_inp(4,1:30000),lineStyleDefault{2},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{4})
xlim([0 60])
ylim([-3 3])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$u_\mathrm{FE}$$ (Nm)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;




