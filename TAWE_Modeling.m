PathSetup;
PathSetup;

%% Declare Time Variable and System
t=sym('t','real');
System=kSystem('TAWE',t);
WORLD=System.Space.RootFrame;
System.Space.setPrecompile(true);
System.Model.setPrecompile(true);

%% Declare Constant Parameters
C=System.genParam({
             'empty' 'empty' 0; 
             });
                     
%% Declare Continuous States (Generalized Coordinates) 
[q1dev,q1flex,q1sup]...
=System.genCont({
            'q1dev' 'Forearm Deviation Angle';
            'q1flex' 'Forearm Flexion Angle';
            'q1sup' 'Forearm Supination Angle';
            });
        
[q2w1,q2w2,q2w3,q2w4,q2w5,q2w6]...
=System.genCont({
            'q2w1' 'Exoskeleton Wrist Link Angle 1';
            'q2w2' 'Exoskeleton Wrist Link Angle 2';
            'q2w3' 'Exoskeleton Wrist Link Angle 3';
            'q2w4' 'Exoskeleton Wrist Link Angle 4';
            'q2w5' 'Exoskeleton Wrist Link Angle 5';
            'q2w6' 'Exoskeleton Wrist Link Angle 6';
            });
        
%% Declare System Discrete Variables 
P_Sys=System.genDisc({
             'g' 'Gravity Acceleration';
             'gUncertain' 'Uncertain Body Gravity Acceleration';
             'c' 'Global Soft Constraint Factor';
             'b1act' 'Arm Damping Factor';
             'b2act' 'Exoskeleton Damping Factor';
             'rTetra' 'Simplex3 Radius'
             });
         
P_Arm=System.genDisc({
             'l11x' 'Joint Deviation along x from IMU1 axis (after rotation)';
             'l11y' 'Joint Deviation along y from IMU1 axis (after rotation)';
             'l11z' 'Joint Deviation along z from IMU1 axis (after rotation)';
             'l12x' 'Hand IMU Deviation along x from flex frame (after rotation)';
             'l12y' 'Hand IMU Deviation along y from flex frame (after rotation)';
             'l12z' 'Hand IMU Deviation along z from flex frame (after rotation)';
             'l13x' 'Connection Point Deviation along x from IMU2 frame (after rotation)';
             'l13y' 'Connection Point Deviation along y from IMU2 frame (after rotation)';
             'l13z' 'Connection Point Deviation along z from IMU2 frame (after rotation)';
             });
         
P_Exo=System.genDisc({
             'l21x' 'Exo Device Distance along x from IMU1 axis';
             'l21y' 'Exo Device Distance along y from IMU1 axis';
             'l21z' 'Exo Device Distance along z from IMU1 axis';   
             'l21' 'Linkage 1 Length';
             'l22' 'Linkage 2 Length';
             'l23' 'Linkage 3 Length';
             'l24' 'Linkage 4 Length';
             'l25' 'Linkage 5 Length';
             });
         
P_Init1=System.genDisc({
             'rf1IMU1x' 'Initial IMU1 Orientation x';
             'rf1IMU1y' 'Initial IMU1 Orientation y';
             'rf1IMU1z' 'Initial IMU1 Orientation z';
             'rf1Devx' 'Initial Wrist Base Orientation x with respect to rotated supination axis';
             'rf1Devy' 'Initial Wrist Base Orientation y with respect to rotated supination axis';
             'rf1Devz' 'Initial Wrist Base Orientation z with respect to rotated supination axis';
             'rf1IMU2x' 'Initial IMU2 Orientation x';
             'rf1IMU2y' 'Initial IMU2 Orientation y';
             'rf1IMU2z' 'Initial IMU2 Orientation z';
             });
         
P_Init2=System.genDisc({
             'rf2Medx' 'Initial Device Orientation x (Forearm)';
             'rf2Medy' 'Initial Device Orientation y (Forearm)';
             'rf2Medz' 'Initial Device Orientation z (Forearm)';
             'riw1' 'w1 Initial Displacement';
             'riw2' 'w2 Initial Displacement';
             'riw3' 'w3 Initial Displacement';
             'riw4' 'w4 Initial Displacement';
             'riw5' 'w5 Initial Displacement';
             'riw6' 'w6 Initial Displacement';
             });
         
P_COM1=System.genDisc({
             'lm13x' 'Hand COM x';
             'lm13y' 'Hand COM y';
             'lm13z' 'Hand COM z';
             });
         
P_COM2=System.genDisc({
             'lm22x' 'Exoskeleton Pan COM x';
             'lm22y' 'Exoskeleton Pan COM y';
             'lm22z' 'Exoskeleton Pan COM z';
             'lm23x' 'Exoskeleton Link 1 COM x';
             'lm23y' 'Exoskeleton Link 1 COM y';
             'lm23z' 'Exoskeleton Link 1 COM z';
             'lm24x' 'Exoskeleton Link 2 COM x';
             'lm24y' 'Exoskeleton Link 2 COM y';
             'lm24z' 'Exoskeleton Link 2 COM z';
             'lm25x' 'Exoskeleton Link 3 COM x';
             'lm25y' 'Exoskeleton Link 3 COM y';
             'lm25z' 'Exoskeleton Link 3 COM z';
             'lm26x' 'Exoskeleton Link 4 COM x';
             'lm26y' 'Exoskeleton Link 4 COM y';
             'lm26z' 'Exoskeleton Link 4 COM z';
             'lm27x' 'Exoskeleton Link Joint COM x';
             'lm27y' 'Exoskeleton Link Joint COM y';
             'lm27z' 'Exoskeleton Link Joint COM z';
             });
         
P_Inert1=System.genDisc({
             'm13' 'Hand  Mass';
             'i13xx' 'Hand Inertia xx';
             'i13xy' 'Hand Inertia xy';
             'i13xz' 'Hand Inertia xz';
             'i13yy' 'Hand Inertia yy';
             'i13yz' 'Hand Inertia yz';
             'i13zz' 'Hand Inertia zz';
             });
        
P_Inert2=System.genDisc({
             'm22' 'Exoskeleton Pan Mass';
             'i22xx' 'Exoskeleton Pan Inertia xx';
             'i22xy' 'Exoskeleton Pan Inertia xy';
             'i22xz' 'Exoskeleton Pan Inertia xz';
             'i22yy' 'Exoskeleton Pan Inertia yy';
             'i22yz' 'Exoskeleton Pan Inertia yz';
             'i22zz' 'Exoskeleton Pan Inertia zz';
             'm23' 'Exoskeleton Link 1 Mass';
             'i23xx' 'Exoskeleton Link 1 Inertia xx';
             'i23xy' 'Exoskeleton Link 1 Inertia xy';
             'i23xz' 'Exoskeleton Link 1 Inertia xz';
             'i23yy' 'Exoskeleton Link 1 Inertia yy';
             'i23yz' 'Exoskeleton Link 1 Inertia yz';
             'i23zz' 'Exoskeleton Link 1 Inertia zz';
             'm24' 'Exoskeleton Link 2 Mass';
             'i24xx' 'Exoskeleton Link 2 Inertia xx';
             'i24xy' 'Exoskeleton Link 2 Inertia xy';
             'i24xz' 'Exoskeleton Link 2 Inertia xz';
             'i24yy' 'Exoskeleton Link 2 Inertia yy';
             'i24yz' 'Exoskeleton Link 2 Inertia yz';
             'i24zz' 'Exoskeleton Link 2 Inertia zz';
             'm25' 'Exoskeleton Link 3 Mass';
             'i25xx' 'Exoskeleton Link 3 Inertia xx';
             'i25xy' 'Exoskeleton Link 3 Inertia xy';
             'i25xz' 'Exoskeleton Link 3 Inertia xz';
             'i25yy' 'Exoskeleton Link 3 Inertia yy';
             'i25yz' 'Exoskeleton Link 3 Inertia yz';
             'i25zz' 'Exoskeleton Link 3 Inertia zz';
             'm26' 'Exoskeleton Link 4 Mass';
             'i26xx' 'Exoskeleton Link 4 Inertia xx';
             'i26xy' 'Exoskeleton Link 4 Inertia xy';
             'i26xz' 'Exoskeleton Link 4 Inertia xz';
             'i26yy' 'Exoskeleton Link 4 Inertia yy';
             'i26yz' 'Exoskeleton Link 4 Inertia yz';
             'i26zz' 'Exoskeleton Link 4 Inertia zz';
             'm27' 'Exoskeleton Link Joint Mass';
             'i27xx' 'Exoskeleton Link Joint Inertia xx';
             'i27xy' 'Exoskeleton Link Joint Inertia xy';
             'i27xz' 'Exoskeleton Link Joint Inertia xz';
             'i27yy' 'Exoskeleton Link Joint Inertia yy';
             'i27yz' 'Exoskeleton Link Joint Inertia yz';
             'i27zz' 'Exoskeleton Link Joint Inertia zz';
             });
         
P_Uncertain=System.genDisc({
                             'm12' 'Hand Load Mass';
                             'i12xx' 'Hand Load Inertia xx';
                             'i12xy' 'Hand Load Inertia xy';
                             'i12xz' 'Hand Load Inertia xz';
                             'i12yy' 'Hand Load Inertia yy';
                             'i12yz' 'Hand Load Inertia yz';
                             'i12zz' 'Hand Load Inertia zz';
                             });

U1=System.genInput({
             'u1dev' 'Arm Deviation Input';
             'u1flex' 'Arm Flexion Input';
             });

U2=System.genInput({
             'u21' 'Exo Skeleton Input 1';
             'u22' 'Exo Skeleton Input 2';
             });
         
%These are actually parameters. We use the toolbox to generate the numerical Jacobians of these parameters.
U_Uncertain=System.genInput({
                             'handmpt1' 'Hand Load Uncertain Mass 1'
                             'handmpt2' 'Hand Load Uncertain Mass 2'
                             'handmpt3' 'Hand Load Uncertain Mass 3'
                             'handmpt4' 'Hand Load Uncertain Mass 4'
                             }); 
         
         
[xx3,yy3,zz3]=System.genNHSignal({
                                 'xx3' 'Euler3 X' 1;
                                 'yy3' 'Euler3 Y' 1;
                                 'zz3' 'Euler3 Z' 1;
                                });
        
%% Kinematic and Dynamic Modeling
TAWE_Modeling_Kinematics;
TAWE_Modeling_Dynamics;

%% Compile Dynamics
System.Model.Init('Full');
[FixArm]=System.Model.genNode({'FixArm'});
System.Model.plotGraph(2);
FixArm.genFlowEOM({...
                   'AttachExoEndDx','AttachExoEndDy','AttachExoEndDz',...
                   'AttachExoEndRx','AttachExoEndRy','AttachExoEndRz',...
                   'SupConsY',...
                   }.');
System.Model.makeDynamics('num');