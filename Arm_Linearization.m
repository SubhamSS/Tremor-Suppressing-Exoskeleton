%% Linearization of the Arm model

% Simulation code for "Model Predictive Control for Tremor Suppressing Exoskeleton"
% by Subham Samal (subhamsamal@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech

% The authors thank Jiamin Wang (jmechw@vt.edu) for their constant guidance
% and support throughout the project

% This script contains the formulation to obtain the linearized model of
% the non-linear Arm Dynamical System. The code generates mex files 
% which are then used during simulation 

% mkdir Linearization_ARM
PathSetup;
%% Load the dynamical model

load ./Dynamics_Arm/BaseEOM_Arm.mat
load ./Dynamics_Arm/ModelInfo_Arm.mat

%% Goal of linearization
% xxdot = ff(xx) = ff(xx_0) + AA(xx_0)*(xx-xx_0)
% leading to the discrete time form
% xx_t+i+1 = xx_t+i + hh*ff(xx_t) + hh*AA(xx_t)*(xx_t+i - xx_t) 
% where AA is the state jacobian
% xx_t is the current state, i starts from 0

% Hence, solving AA will lead to efficient approximation for MPC
% calculation
%%
Arm_Parameters;
%%  Load Matrices
disp("Preparing Matrices...")

syms t positive;

MM_lin = BaseEOM.InertialMatrix; % Inertia Matrix
CC_lin = BaseEOM.CoriolisMatrix+BaseEOM.CentripetalMatrix; % Coriolis and Centripetal Matrix
hh_lin = BaseEOM.CoriolisForce+...% Generalized Forces
    BaseEOM.DampingForce+...
    BaseEOM.ExternalForce+...
    BaseEOM.InputForce+...
    BaseEOM.PotentialForce;

JuT_lin = ModelInfo.Dynamics.InputJacobian; % Input Jacobian
uu_lin = ModelInfo.Variables.Input; % Control Input
disp("done!")

%% Substitute Model Parameters

disp("Substituting Model Parameters...")
pp = symvar(hh_lin).';
pp = [pp(1:3, :);pp(8:27, :);pp(34:43, :)];
ppVal = [ppValCtrlNoWKI(4);ppValCtrlNoWKI(1:2);ppValCtrlNoWKI(36:41);ppValCtrlNoWKI(29:34);...
    ppValCtrlNoWKI(10:12);ppValCtrlNoWKI(25:27);ppValCtrlNoWKI(35);ppValCtrlNoWKI(28);ppValCtrlNoWKI(6);...
    ppValCtrlNoWKI(19:21);ppValCtrlNoWKI(16:18);ppValCtrlNoWKI(22:24)];

MM_lin1 = vpa(subs(MM_lin,[pp],[ppVal]));
CC_lin1 = vpa(subs(CC_lin,[pp],[ppVal]));
hh_lin1 = vpa(subs(hh_lin,[pp],[ppVal]));
JuT_lin1 = vpa(subs(JuT_lin,[pp],[ppVal]));
disp("done!")

%% Constrained System 
% Define States 
disp("Preparing Constrained System...")

qqdev_ = sym("q1dev__dt_0_"); % Wrist Radial Ulnar Deviation (RUD) Rotation on z axis in the wrist frame
qqflex_ = sym("q1flex__dt_0_"); % Wrist Flexion Extension (FE) Rotation on x axis in the wrist frame
qqsup_ = sym("q1sup__dt_0_"); % Wrist Pronation Supination (SP), which is constrained to FE and
                              % RUD, Rotation on y axis in the wrist frame
qqf = [qqdev_;qqflex_];

% Angular Velocities
qqfdot = pDiff(qqf, t);

% State Vector (unconstrained)
xx = [qqf;qqfdot];

% All States Vector
states = [qqdev_;qqflex_;qqsup_;pDiff(qqdev_,t);pDiff(qqflex_,t);pDiff(qqsup_,t)];

% Inputs
inputs = uu_lin;
uu1lin = uu_lin(1:2);

% M_dot
MM_lin1dot = pDiff(MM_lin1,t);

% JlamMapFull and JlamMapFulldot placeholder matrices
JL = sym("jl%d",[3,2]);
JLdot = sym("jld",[3,2]);

% Calculated properties of constrained system
MMRC = JL.'*MM_lin1;
MMC = MMRC*JL;

MMdot = JL.'*MM_lin1dot*JL;

CCC = JL.'*(CC_lin1)*JL ...
                    + MMRC*JLdot; 
hhC = JL.'*hh_lin1; 
JacU1Tpose = JL.'*JuT_lin1(:,1:2); % Jacobian of control input (transposed)

% Model equation
ZZC = (-CCC*qqfdot-hhC+JacU1Tpose*uu1lin);

%% Numerical-Efficient Solution

% To substitute M_dot numerically
WW = sym("ww%d",[2,2]);
% To substitute ZZC numerically
zz = sym("zz%d",[2,1]);

% AA_21 (AA_21 is a 2*2 matrix calculated as AA21_a +AA21_b
AA21_a = jacobian(-WW*MMdot*WW*zz,qqfdot);
disp("Generated AA21_a")

AA21b_11 = pDiff(ZZC(1),qqf(1));
AA21b_12 = pDiff(ZZC(1),qqf(2));
AA21b_21 = pDiff(ZZC(2),qqf(1));
AA21b_22 = pDiff(ZZC(2),qqf(2));
disp("Generated AA21_b")

AA22_11 = pDiff(ZZC(1),qqfdot(1));
AA22_12 = pDiff(ZZC(1),qqfdot(2));
AA22_21 = pDiff(ZZC(2),qqfdot(1));
AA22_22 = pDiff(ZZC(2),qqfdot(2));
disp("Generated AA22")
%% Generate function script
% This could take some time         

AA21_afun = matlabFunction(AA21_a, 'vars', {states,JLdot, JL, WW,zz},"File","./Linearization_ARM/AA21_a_Val_Arm");
disp("Generated matlabfunction for AA21_a")

AA21b11fun = matlabFunction(AA21b_11, 'vars', {states,inputs,JLdot, JL},"File","./Linearization_ARM/AA21b_11_Val_Arm");
disp("Generated matlabfunction for AA21b_11")
AA21b12fun = matlabFunction(AA21b_12, 'vars', {states,inputs,JLdot, JL},"File","./Linearization_ARM/AA21b_12_Val_Arm");
disp("Generated matlabfunction for AA21b_12")
AA21b21fun = matlabFunction(AA21b_21, 'vars', {states,inputs,JLdot, JL},"File","./Linearization_ARM/AA21b_21_Val_Arm");
disp("Generated matlabfunction for AA21b_21")
AA21b22fun = matlabFunction(AA21b_22, 'vars', {states,inputs,JLdot, JL},"File","./Linearization_ARM/AA21b_22_Val_Arm");
disp("Generated matlabfunction for AA21b_22")

AA2211fun = matlabFunction(AA22_11, 'vars', {states,JLdot, JL},"File","./Linearization_ARM/AA22_11_Val_Arm");
disp("Generated matlabfunction for AA22_11")
AA2212fun = matlabFunction(AA22_12, 'vars', {states,JLdot, JL},"File","./Linearization_ARM/AA22_12_Val_Arm");
disp("Generated matlabfunction for AA22_12")
AA2221fun = matlabFunction(AA22_21, 'vars', {states,JLdot, JL},"File","./Linearization_ARM/AA22_21_Val_Arm");
disp("Generated matlabfunction for AA22_21")
AA2222fun = matlabFunction(AA22_22, 'vars', {states,JLdot, JL},"File","./Linearization_ARM/AA22_22_Val_Arm");
disp("Generated matlabfunction for AA22_22")

%% Generate precompiled functions 
codegen ./Linearization_ARM/AA21_a_Val_Arm -args {zeros(18,1),zeros(9,2),zeros(9,2),zeros(2,2),zeros(2,1)}  -o ./Linearization_ARM/AA21_a_Val_Arm_mex -report

codegen ./Linearization_ARM/AA21b_11_Val_Arm -args {zeros(18,1),zeros(8,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA21b_11_Val_Arm_mex -report
codegen ./Linearization_ARM/AA21b_12_Val_Arm -args {zeros(18,1),zeros(8,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA21b_12_Val_Arm_mex -report
codegen ./Linearization_ARM/AA21b_21_Val_Arm -args {zeros(18,1),zeros(8,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA21b_21_Val_Arm_mex -report
codegen ./Linearization_ARM/AA21b_22_Val_Arm -args {zeros(18,1),zeros(8,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA21b_22_Val_Arm_mex -report

codegen ./Linearization_ARM/AA22_11_Val_Arm -args {zeros(18,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA22_11_Val_Arm_mex -report
codegen ./Linearization_ARM/AA22_12_Val_Arm -args {zeros(18,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA22_12_Val_Arm_mex -report
codegen ./Linearization_ARM/AA22_21_Val_Arm -args {zeros(18,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA22_21_Val_Arm_mex -report
codegen ./Linearization_ARM/AA22_22_Val_Arm -args {zeros(18,1),zeros(9,2),zeros(9,2)}  -o ./Linearization_ARM/AA22_22_Val_Arm_mex -report

