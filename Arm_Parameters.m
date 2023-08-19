%% Parameters for the Arm simulations

% Simulation code for "Model Predictive Control for Tremor Suppressing Exoskeleton"
% by Subham Samal (subhamsamal@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech

% The authors thank Jiamin Wang (jmechw@vt.edu) for their constant guidance
% and support throughout the project

% This script contains the simulation, model, and 
% control parameters for all Arm simulations. 

%% Simulation Parameters
simFreq=500; % Simulation sample rate

hh=1/simFreq; % Simulation step size
tt=0; % Time variables
flowNum=1; % This is used to select different sets of model constraints, set to 1 for all TAWE simulations
tSpan=0:hh:150; % Default time span

%% Model Parameters

% System parameters
rTetra=0.2/3*sqrt(6); % Dimension of the load tetrahedron (for center of mass estimation) on the hand

% The parameters of P_Sys are: 
% (1) Gravity acceleration
% (2) Gravity acceleration of the uncertain body (set to 1 for simulation, 
%     set to 0 when load tetrahedron (for center of mass estimation) is used)
% (3) Soft constraint Stiffness (used to prevent accumulated drifting error
%     during numerical integration of constrained system. This term will
%     not be active if the constraint is not significantly violated)
% (4) Hand damping
% (5) Exoskeleton damping
% (6) Load tetrahedron (for center of mass estimation) dimension parameter.

P_Sys_Val=  [9.8067;  9.8067;  100;  0.1;  0.05;  rTetra]; 

% Kinematic parameters obtained from CAD designs, which matches with the 3D
% models in visualization. Modifications are not recommended.

% Constant angles in the forearm/wrist model
P_Init1_Val=1*pi*[0;0;0; -10;15;5; 0;0;-10]/180; 

% Translational displacements in the wrist model obtained from CAD
P_Arm_Val=[-6.28; 70; -41; ...
           11.76; -60.87; -27.98;...
           -5; -30; 1]*1e-3; 
P_Arm_Val(4:6)=-ypr2Mat(P_Init1_Val(7:9))*P_Arm_Val(4:6); % A parameter adjustment to match with 3D model

% Hand Center of Mass displcaement parameters obtained from CAD
P_COM1_Val=[
            -15.38; 16.49; -28.05
            ]*1e-3;
        
% Hand mass and moment of inertia obtained from CAD
% columns: (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
P_Inert1=[
          0.656;    0.656*[+2.259e-03;-5.823e-04;+1.027e-04;+1.430e-03;+2.206e-04;+3.354e-03];
          ];

% Uncertain load mass and inertia (set to zeros for now, to be used in future works)
P_Uncertain_Val=zeros(7,1);
P_Uncertain_Ctrl=zeros(7,1);
      
% All parameters used for simulation, controller without WKI, and
% controller with WKI
ppValSim=[P_Sys_Val;P_Arm_Val;P_Init1_Val;P_COM1_Val;P_Inert1;P_Uncertain_Val];
ppValCtrlNoWKI=[P_Sys_Val;P_Arm_Val;P_Init1_Val;P_COM1_Val;P_Inert1;P_Uncertain_Ctrl]; 

ppNum=numel(ppValSim); % Total model parameter count


