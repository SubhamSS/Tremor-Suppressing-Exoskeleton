%% The Wearable Wrist Exoskeleton (TAWE) Control Simulation with tremor suppression (Figs 7,8)

% Simulation code for "Model Predictive Control for Tremor Suppressing Exoskeleton"
% by Subham Samal (subhamsamal@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech

% The authors thank Jiamin Wang (jmechw@vt.edu) for their constant guidance
% and support throughout the project

% This script contains the simulation code for the wearable wrist
% exoskeleton for tremor suppression using MPC and BMFLC 

PathSetup;

%% 3D Visualization (Uncomment to load 3D models)
% Below is the initialization codes for the 3D visualization of the
% wrist exoskeleton, which includes 3D models and main coordinate frames. 
% Uncomment the code in this section to initiate visualization. (Requires OpenGL for smooth animation)

Sim=sAxes('TAWE Simulation',3,'Path_TAWE.mat',@numTF_TAWE_mex);
Sim.setAxesProp('ArmIMU1Frame',1*[-0.2 -0.1 -0.2;0.2 0.5 0.2],[210 15]).setPresetTF('WORLD');

[ArmPose,ExoPose,HandCOMTraj,HandTetra]...
    =Sim.genPlot({'ArmPose';'ExoPose';'HandCOMTraj';'HandTetra'});
ArmPose.setLineSpec('-','none','r',2).setPlotPoint({'ArmIMU1Frame';'ArmDevFlexFrame';'ArmIMU2Frame'});
ExoPose.setLineSpec('-','none','b',2).setPlotPoint({'ExoIMU1Frame';'ExoMedFrame';'ExoPanFrame';'ExoLink12Frame';...
                                                    'ExoLink23Frame';'ExoLink34Frame';'ExoLink45Frame';'ExoAttachFrame';'ExoIMU2Frame'});
HandTetra.setLineSpec(':','+','g',3).setPlotPoint({'Handmpt1Frame';'Handmpt2Frame';'Handmpt3Frame';'Handmpt4Frame';'Handmpt2Frame';'Handmpt1Frame';'Handmpt3Frame';'Handmpt1Frame';'Handmpt4Frame';});

[ArmBaseCAD,ArmHandCAD...
     ,ExoBaseCAD,ExoPanCAD,ExoLink1CAD,ExoLink2CAD,...
     ExoLink3CAD,ExoLink4CAD,ExoLinkJointCAD...
    ]...
=Sim.genPatch({
               'ArmBaseCAD'  'ArmIMU1Frame';...
               'ArmHandCAD'  'ArmIMU2AdjustFrame';...
               'ExoBaseLink'  'ExoIMU1Frame';...
               'ExoPan'  'ExoPanFrame';...
               'ExoLink1'  'ExoLink12Frame';...
               'ExoLink2'  'ExoLink23Frame';...
               'ExoLink3'  'ExoLink34Frame';...
               'ExoLink4'  'ExoLink45Frame';...
               'ExoLinkJoint'  'ExoAttachFrame';...
               });

ArmBaseCAD.setFaceProp([0.8 0.5 0.5],0.5).setModel('taweArmBase.STL',0.1,[],0.001);
ArmHandCAD.setFaceProp([0.8 0.5 0.5],0.5).setModel('taweArmHand.STL',0.1,[],0.001);
ExoBaseCAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoBaseLink.STL',0.1,[],1);
ExoPanCAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoPan.STL',0.1,[],1);
ExoLink1CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink1.STL',0.025,[],0.001);
ExoLink2CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink2.STL',0.025,[],0.001);
ExoLink3CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink3.STL',0.025,[],0.001);
ExoLink4CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink4.STL',0.025,[],0.001);
ExoLinkJointCAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLinkJoint.STL',1,[],0.001);
Sim.Patch(1).setFaceProp([0.8 0.8 1],1).setModel([],1,[],0.1);

%% Load Model Parameters (Do not comment)

TAWE_Parameters; % This loads the TAWE simulation parameters
%%

% This loads the tracking reference (position, velocity, acceleration). The
% reference can be generated in TAWE_Simulation_Setup
load('simTAWERef','tremorData');
% For this simulation, tracking reference is kept as stationary (0 in both FE and RUD)
refData = zeros(2,75001);
refdtData = zeros(2,75001);
refddtData = zeros(2,75001);

%% Simulation Initialization

tSpan1 = 0:hh:60; %Simulation Time Span

% State vector qq and qqdot, which are:
% (1) Wrist Radial Ulnar Deviation (RUD) (rad) Rotation on z axis in the wrist frame
% (2) Wrist Flexion Extension (FE) (rad) Rotation on x axis in the wrist frame
% (3) Wrist Pronation Supination (SP) (rad), which is constrained to FE and
%     RUD, Rotation on y axis in the wrist frame
% (4-9) TAWE Mechanism Joint Angles (rad)
simStateVal=zeros(18,1);
simParamVal=ppValSim; % Simulation  parameters
ctrlParamVal=ppValCtrlNoWKI; % Controller Parameters

% The input vectors uu are defined as 
% (1-2) human inputs at the wrist FE and RUD directions
% (3-4) control input from TAWE at the first two links
% (5-8) virtual input at the load tetrahedron vertices (-z direction) 
%       (for center of mass estimation), used to acquire the Jacobian for 
%       uncertain load parameter more conveniently. Considered 0 for this
%       simulation (to be explored in future works)

simInputVal=zeros(8,1); % Used for simulation 
ctrllnputVal=zeros(8,1); % Used for controller calculation only
taweInput=zeros(2,1); % The TAWE input term

consForce=zeros(7,1); % Constraint force for simulation process (not used in controller design)

% Nonholonomic state used in constraint only. Note that are not internal 
% states of the system, but Euler Angles qq_nh that satisfies
% RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
% Hence, the rotation by qq1 is split to two identical rotations by qq_nh
% in sequance. The wrist kinematic model uses this to form a z-x-z-x 
% sequenced rotation joint, where all the z angles (on RUD direction) are the same, and all
% the x angles (on FE direction) are the same. This assumes that the wrist rotation is
% devided evenly into the radiocarpal and midicarpal joints.
nhSignal=zeros(3,1);

simStateData=zeros(18,numel(tSpan1)); % Store state data from simulation
simInputData=zeros(8,numel(tSpan1)); %store inputs

%% Simulation Loop
tic;
for ii=1:numel(tSpan1)

    % just for checkpoint
    if rem(ii,100)==1
        disp(ii+" of "+numel(tSpan1)+" steps completed")
    end
    
    tt=tSpan1(ii); % Update time variable
    simStateData(:,ii)=simStateVal(1:end); % Store qq
    qq=simStateVal(1:end/2);
    qqdot=simStateVal(end/2+1:end);
    
    ref=refData(:,ii);
    refdt=refdtData(:,ii);
    refddt=refddtData(:,ii);
        
    % Get System Properties
    [~,~,~,~,~,Quat,MM,~,GFF,~,JacU,JacCons1,~,~,CenMat,InertLeftComponent,InertRightComponent1]...
        =System_TAWE_mex(tSpan1(ii),simStateVal,ctrlParamVal,ctrllnputVal,nhSignal);
    MM=sum(MM,3); % Sum up inertia matrix of all bodies
    JacU=sum(JacU,3).'; % Sum up constraint Jacobian Matrices of all inputs
    JacCons1=sum(JacCons1,3); % Sum up constraint Jacobian Matrices of all constraints
    GFF=sum(GFF,2); % Sum up generalized forces 
    CenMat=sum(CenMat,3); % Sum up all centripetal matrices'

    % Properties used to numerically calculate Coriolis matrix and
    % constraint Jacobian derivative
    simStateVal2=simStateVal;
    % Use the current state velocity to estimate the position after a
    % short time step
    simStateVal2(1:end/2)=simStateVal2(1:end/2)+simStateVal(end/2+1:end)*(hh);

    % Calculate the equivalent qq_nh used in the constraint (see
    % above), Note that are not internal states of the system, but Euler
    % Angles qq_nh that satisfies
    % RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
    % Hence, the rotation by qq1 is split to two identical rotations by qq_nh
    % in sequance.
    rotQuat=ypr2Quat(simStateVal2([2 3 1])); % Unit quaternion equivalent to wrist rotation
    halfCos=sqrt((rotQuat(1)+1)/2);
    rotQuatSquareRoot=[halfCos;rotQuat(2:4)/2/halfCos]; % Half of the rotQuat rotation
    nhSignal2=quat2YPR(rotQuatSquareRoot)*2; % Calculate qq_nh estimate after a short time step
    % Obtain estimated properties after a short time step
    [~,~,~,~,~,~,~,~,~,~,~,JacCons2,~,~,~,~,InertRightComponent2]=...
        System_TAWE_mex(tSpan1(ii),simStateVal2,ctrlParamVal,ctrllnputVal,nhSignal2);

    % Formulate constraint properties (please also see supplemental document)
    Jlam1=JacCons1(:,1:2); % The Jacobian for preseved coordinates of the constrained model
    Jlam2=JacCons1(:,3:end); % The Jacobian for internal coordinates of the constrained model
    JlamMap=-Jlam2\Jlam1; % The mapper between preseved coordinates and internal coordinates
    JlamMapFull1=[eye(2);JlamMap]; % The mapper between preseved coordinates and all coordinates

    % This step is similar to above, except using the estimated
    % properties
    Jlam1=JacCons2(:,1:2);
    Jlam2=JacCons2(:,3:end);
    JlamMap=-Jlam2\Jlam1;
    JlamMapFull2=[eye(2);JlamMap];

    % Numerically calculate derivatives of matrices (the Jacobian
    % matrices do not involve velocity (i.e. qqdot))
    JlamMapFull=JlamMapFull1;
    JlamMapFullDot=(JlamMapFull2-JlamMapFull1)/(hh);
    InertRightComponent=InertRightComponent1;
    InertRightComponentDot=(InertRightComponent2-InertRightComponent1)/(hh);

    % Numerically calculate coriolis matrix
    CorMat=InertLeftComponent*InertRightComponentDot;

    % Calculated properties of constrained system
    MMRightCombine = JlamMapFull.'*MM;
    MMCombine = MMRightCombine*JlamMapFull; % Inertia matrix
    CCCombine = JlamMapFull.'*(CorMat+CenMat)*JlamMapFull ...
        + MMRightCombine*JlamMapFullDot; % Coriolis and Centripetal Matrix
    GFFCombine = JlamMapFull.'*GFF; % Generalized force
    JacU1Tpose = JlamMapFull.'*JacU(1:2,:).'; % Jacobian of user input (transposed)
    JacU2Tpose = JlamMapFull.'*JacU(3:4,:).'; % Jacobian of TAWE input (transposed)

    % values for linearization
    statesval = [simStateVal];
    inputsval = [0;0;simInputVal(3:8)];
    JLVal = JlamMapFull;
    JLValdot = JlamMapFullDot;
    invMM = inv(MMCombine);
    ZZVal = -CCCombine*simStateVal(10:11)-(GFFCombine)+JacU2Tpose*simInputVal(3:4);
    
    % obtain state space matrices using linearization
    step = 1;
    if rem(ii,5)==1 % Linearizing every 5th step
        AA_21a = AA21_a_Val_mex(statesval,JLValdot,JLVal,invMM,ZZVal);
        AA_21b = invMM*[AA21b_11_Val_mex(statesval,inputsval,JLValdot,JLVal),AA21b_12_Val_mex(statesval,inputsval,JLValdot,JLVal);...
            AA21b_21_Val_mex(statesval,inputsval,JLValdot,JLVal),AA21b_22_Val_mex(statesval,inputsval,JLValdot,JLVal)];
        AA_21 = AA_21a + AA_21b;
        AA_22 = invMM*[AA22_11_Val_mex(statesval,JLValdot,JLVal),AA22_12_Val_mex(statesval,JLValdot,JLVal);...
            AA22_21_Val_mex(statesval,JLValdot,JLVal),AA22_22_Val_mex(statesval,JLValdot,JLVal)];
        Ac = [zeros(2,2), eye(2,2); AA_21, AA_22];
        Ad = eye(4)+ hh*Ac*step;
        Bc = [zeros(2,2);invMM*JacU2Tpose];
        Bd = Bc*hh*step;
    end
    Consterm = hh*step*([simStateVal(10:11);invMM*ZZVal]-Ac*[simStateVal(1:2);simStateVal(10:11)]-Bc*simInputVal(3:4));
    X_in = [simStateVal(1:2);simStateVal(10:11)]; 
    
    hor=20;            % Horizon
    path2 = [refData(:,ii:step:ii+hor*step-1);refdtData(:,ii:step:ii+hor*step-1)]; %reference path

    % State error weight
    Q=diag([200 1000 0.1 0.1]);
    % final state error weight
    P = Q;     
    % input weights
    R=0.4*eye(2); 

    % Input constraints 
    uMin=[-9.5;-12.2];      
    uMax=[11;7.1];
    
    % Ang pos and vel constraints
    xx_lb = [-0.8;-1.5;-3*pi/2;-3*pi/2]; 
    xx_ub = [0.45;1.5;3*pi/2;3*pi/2];

    %Obtain optimal input from MPC Controller
    taweInput = mpc_no_BMFLC(Ad,Bd,P,Q,R,hor,X_in,Consterm,path2, uMax, uMin, xx_lb, xx_ub, hh);

    simInputVal(3:4)=taweInput(1:2); % Insert MPC input for simulation

    % Tremor excitation from user
    userInput=tremorData(:,ii); 
    simInputVal(1:2)=userInput;

    %Store input data
    simInputData(:,ii)=simInputVal(1:end);
    
    % Using Runge Kutta 4th method to integrate system ODE
    [simStateVal,nhSignal,consForce]=rk4Hybrid(@Flow_TAWE_mex,hh,flowNum,tSpan1(ii),...
                                        simStateVal,simParamVal,simInputVal,nhSignal,consForce);

    
    % The codes below refreshes the 3D visualization for run-time animation
    % Uncomment to run.
    if rem(ii,15)==0 %adjust the frame rate through the number
        Sim.drawNow(tSpan1(ii),simStateVal,simParamVal,simInputVal,nhSignal); % requires 3D environment
    end
end
toc;
%% Store Simulation data to plot (Uncomment to Run)
 
save(strcat('./Data/zero_track_TAWE_MPC'),'simStateData') % Save result
save(strcat('./Data/zero_track_TAWE_MPC_inp'),'simInputData')

