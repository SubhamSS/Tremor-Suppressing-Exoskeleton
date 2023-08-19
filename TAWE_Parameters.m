%% Parameters for the Wearable Wrist Exoskeleton (TAWE) simulations

% Simulation code for "Model Predictive Control for Tremor Suppressing Exoskeleton"
% by Subham Samal (subhamsamal@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech

% The authors thank Jiamin Wang (jmechw@vt.edu) for their constant guidance
% and support throughout the project

% This script contains the simulation, model, and control parameters
% for all wearable wrist exoskeleton (TAWE) simulations. 

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

% Constant angles/offsets in TAWE
P_Init2_Val=[0;0;0;-0.3699;0.7732;-1.8355;0.7874;0.1527;0.3268];

% Translational displacements in the wrist model obtained from CAD
P_Arm_Val=[-6.28; 70; -41; ...
           11.76; -60.87; -27.98;...
           -5; -30; 1]*1e-3; 
P_Arm_Val(4:6)=-ypr2Mat(P_Init1_Val(7:9))*P_Arm_Val(4:6); % A parameter adjustment to match with 3D model

% TAWE mechanism link parameters obtained from CAD (known)
P_Exo_Val=[-68.6; 12; 21.5; ...
           120;120;30;10;-25]*1e-3;
P_Exo_Ctrl=P_Exo_Val;

% Hand Center of Mass displcaement parameters obtained from CAD
P_COM1_Val=[
            -15.38; 16.49; -28.05
            ]*1e-3;
        
% TAWE mechanism Center of Mass displcaement parameters obtained from CAD (known)
P_COM2_Val=[
             -3.08; -10.54; -1.87;
            20.42; 48.05; 0;
            34.56; 60; 0;
            9.06; 8.28; 0;
            1.72; -4.49; -10.65;
            -1.32; -18.35; -2.32;
            ]*1e-3;
        
% Hand mass and moment of inertia obtained from CAD
% columns: (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
P_Inert1=[
          0.656;    0.656*[+2.259e-03;-5.823e-04;+1.027e-04;+1.430e-03;+2.206e-04;+3.354e-03];
          ];
      
% TAWE mechanism link masses and moments of inertia obtained from CAD (known)
% columns: (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
P_Inert2=[
          0.092;    0.092*[+2.267e-04;-1.115e-05;-3.116e-06;+1.831e-04;-9.881e-06;+2.478e-04];  ...
          0.028;    0.028*[+1.490e-03;-1.936e-04;-2.869e-04;+1.114e-03;-6.778e-04;+5.672e-04];  ...
          0.032;    0.032*[+2.169e-03;-1.516e-04;+3.318e-04;+1.860e-03;+8.138e-04;+4.813e-04];  ...
          0.012; 	0.012*[+5.682e-05;+2.011e-05;+9.749e-09;+1.630e-04;-2.356e-09;+1.776e-04];  ...
          0.018;    0.018*[+1.602e-04;+4.302e-07;-4.109e-05;+1.493e-04;-2.381e-05;+1.656e-04];  ...
          0.030;    0.030*[+2.260e-04;+1.539e-05;+1.650e-12;+3.106e-04;-6.426e-12;+5.167e-04];  ...
         ];

% Uncertain load mass and inertia (set to zeros for now, to be used in future works)
P_Uncertain_Val=zeros(7,1);
P_Uncertain_Ctrl=zeros(7,1);
      
% All parameters used for simulation and controller 
ppValSim=[P_Sys_Val;P_Arm_Val;P_Exo_Val;P_Init1_Val;P_Init2_Val;P_COM1_Val;P_COM2_Val;P_Inert1;P_Inert2;P_Uncertain_Val];
ppValCtrlNoWKI=[P_Sys_Val;P_Arm_Val;P_Exo_Val;P_Init1_Val;P_Init2_Val;P_COM1_Val;P_COM2_Val;P_Inert1;P_Inert2;P_Uncertain_Ctrl]; 

ppNum=numel(ppValSim); % Total model parameter count

%% Control Parameters

% BMFLC reference model parameters
bmflcRes=16; % Resolution, i.e., frequency components
bmflcRange=(linspace(2,7,bmflcRes)*2*pi).'; % Frequencies from Bandwidth (2-7 Hz by default)
