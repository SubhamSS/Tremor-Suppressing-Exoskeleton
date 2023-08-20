%% Declare Coordinate Frames in Forearm
[ArmIMU1Frame,ArmBaseDevFrame,ArmDevFlexFrame,ArmIMU2Frame,ArmIMU2AdjustFrame....
]...
=System.Space.genNode({
                     'ArmIMU1Frame'
                     'ArmBaseDevFrame';
                     'ArmDevFlexFrame';
                     'ArmIMU2Frame';
                     'ArmIMU2AdjustFrame'
                     });
                 
                 
[ArmHandCOM]...
=System.Space.genNode({
                     'ArmHandCOM';
                     });
                 
                 
                 
[Handmpt1Frame,Handmpt2Frame,Handmpt3Frame,Handmpt4Frame]...
=System.Space.genNode({
                     'Handmpt1Frame';
                     'Handmpt2Frame';
                     'Handmpt3Frame';
                     'Handmpt4Frame';
                     });
                 

%% Declare Links in Forearm
World2AI1F=System.Space.genEdge({'World2AI1F' 'WORLD' 'ArmIMU1Frame'});

q1IMU1Quat=ypr2Quat([P_Init1.rf1IMU1x;P_Init1.rf1IMU1y;P_Init1.rf1IMU1z]);
World2AI1F.setQuat(q1IMU1Quat).genProp;

AI1F2ABDF=System.Space.genEdge({'AI1F2ABDF' 'ArmIMU1Frame' 'ArmBaseDevFrame'});
ABDF2ADFF=System.Space.genEdge({'ABDF2ADFF' 'ArmBaseDevFrame' 'ArmDevFlexFrame'});
ADFF2AI2F=System.Space.genEdge({'ADFF2AI2F' 'ArmDevFlexFrame' 'ArmIMU2Frame'});
AI2F2AI2AF=System.Space.genEdge({'AI2F2AI2AF' 'ArmIMU2Frame' 'ArmIMU2AdjustFrame'});

q1DevAxisQuat=ypr2Quat([P_Init1.rf1Devx;P_Init1.rf1Devy;P_Init1.rf1Devz]);
AI1F2ABDF.setQuat(q1DevAxisQuat).setTransDis([P_Arm.l11x;P_Arm.l11y;P_Arm.l11z]).genProp;

kappa=[q1flex.dot(0);q1sup.dot(0);q1dev.dot(0)];
q1DevFlexQuat=ypr2Quat(kappa);
ABDF2ADFF.setQuat(q1DevFlexQuat).genProp;

ADFF2AI2F.setTransDis([P_Arm.l12x;P_Arm.l12y;P_Arm.l12z]).genProp;

q1IMU2Quat=ypr2Quat([P_Init1.rf1IMU2x;P_Init1.rf1IMU2y;P_Init1.rf1IMU2z]);
AI2F2AI2AF.setQuat(q1IMU2Quat).genProp;


AI2AF2AHCOM=System.Space.genEdge({'AI2AF2AHCOM' 'ArmIMU2AdjustFrame' 'ArmHandCOM'});
AI2AF2AHCOM.setTransDis([P_COM1.lm13x;P_COM1.lm13y;P_COM1.lm13z]).genProp;

Hand2mpt1=System.Space.genEdge({'Hand2mpt1' 'ArmIMU2AdjustFrame' 'Handmpt1Frame'});
Hand2mpt2=System.Space.genEdge({'Hand2mpt2' 'ArmIMU2AdjustFrame' 'Handmpt2Frame'});
Hand2mpt3=System.Space.genEdge({'Hand2mpt3' 'ArmIMU2AdjustFrame' 'Handmpt3Frame'});
Hand2mpt4=System.Space.genEdge({'Hand2mpt4' 'ArmIMU2AdjustFrame' 'Handmpt4Frame'});
Hand2mpt1.setTransDis([-1;-1/sqrt(6);-1/sqrt(3)]*P_Sys.rTetra).genProp;
Hand2mpt2.setTransDis([1;-1/sqrt(6);-1/sqrt(3)]*P_Sys.rTetra).genProp;
Hand2mpt3.setTransDis([0;-1/sqrt(6);2/sqrt(3)]*P_Sys.rTetra).genProp;
Hand2mpt4.setTransDis([0;3/sqrt(6);0]*P_Sys.rTetra).genProp;

%% Declare Coordinate Frames in Exoskeleton
[ExoIMU1Frame,ExoMedFrame,ExoPanFrame,ExoLink12Frame,ExoLink23Frame,ExoLink34Frame,ExoLink45Frame,...
    ExoAttachFrame,ExoIMU2Frame...
]...
=System.Space.genNode({
                     'ExoIMU1Frame'
                     'ExoMedFrame';
                     'ExoPanFrame'
                     'ExoLink12Frame'
                     'ExoLink23Frame'
                     'ExoLink34Frame'
                     'ExoLink45Frame'
                     'ExoAttachFrame'
                     'ExoIMU2Frame'
                     });
                 
                 
[ExoPanCOM,ExoLink1COM,ExoLink2COM,...
    ExoLink3COM,ExoLink4COM,ExoLinkJointCOM]...
=System.Space.genNode({
                     'ExoPanCOM';
                     'ExoLink1COM';
                     'ExoLink2COM';
                     'ExoLink3COM';
                     'ExoLink4COM';
                     'ExoLinkJointCOM';
                     });

                 
%% Declare Links in Exoskeleton
AI1F2EI1F=System.Space.genEdge({'AI1F2EI1F' 'ArmIMU1Frame' 'ExoIMU1Frame'});
AI1F2EI1F.genProp;

EI1F2EMF=System.Space.genEdge({'EI1F2EMF' 'ExoIMU1Frame' 'ExoMedFrame'});
EMF2EPF=System.Space.genEdge({'EMF2EPF' 'ExoMedFrame' 'ExoPanFrame'});
EPF2EL12F=System.Space.genEdge({'EPF2EL12F' 'ExoPanFrame' 'ExoLink12Frame'});
EL12F2EL23F=System.Space.genEdge({'EL12F2EL23F' 'ExoLink12Frame' 'ExoLink23Frame'});
EL23F2EL34F=System.Space.genEdge({'EL23F2EL34F' 'ExoLink23Frame' 'ExoLink34Frame'});
EL34F2EL45F=System.Space.genEdge({'EL34F2EL45F' 'ExoLink34Frame' 'ExoLink45Frame'});
EL45F2EAF=System.Space.genEdge({'EL45F2EAF' 'ExoLink45Frame' 'ExoAttachFrame'});
EAF2EI2F=System.Space.genEdge({'EAF2EI2F' 'ExoAttachFrame' 'ExoIMU2Frame'});

q2MedAxisQuat=ypr2Quat([P_Init2.rf2Medx;P_Init2.rf2Medy;P_Init2.rf2Medz]);
EI1F2EMF.setQuat(q2MedAxisQuat).setTransDis([P_Exo.l21x;P_Exo.l21y;P_Exo.l21z]).genProp;

q2Rot01Quat=ypr2Quat([0;0;P_Init2.riw1+q2w1.dot(0)]);
EMF2EPF.setQuat(q2Rot01Quat).genProp;

q2Rot12Quat=ypr2Quat([P_Init2.riw2+q2w2.dot(0);0;0]);
EPF2EL12F.setQuat(q2Rot12Quat).genProp;

q2Rot23Quat=ypr2Quat([P_Init2.riw3+q2w3.dot(0);0;0]);
EL12F2EL23F.setQuat(q2Rot23Quat).setTransDis([0;P_Exo.l21;0]).genProp;

q2Rot34Quat=ypr2Quat([P_Init2.riw4+q2w4.dot(0);0;0]);
EL23F2EL34F.setQuat(q2Rot34Quat).setTransDis([0;P_Exo.l22;0]).genProp;

q2Rot45Quat=ypr2Quat([0;P_Init2.riw5+q2w5.dot(0);0]);
EL34F2EL45F.setQuat(q2Rot45Quat).setTransDis([0;P_Exo.l23;0]).genProp;

q2Rot56Quat=ypr2Quat([0;0;P_Init2.riw6+q2w6.dot(0)]);
EL45F2EAF.setQuat(q2Rot56Quat).setTransDis([P_Exo.l24;0;P_Exo.l25]).genProp;

q2IMU2Quat=ypr2Quat([P_Init1.rf1IMU2x;P_Init1.rf1IMU2y;P_Init1.rf1IMU2z]);
EAF2EI2F.setQuat(quatConj(q2IMU2Quat)).setTransDis([P_Arm.l13x;P_Arm.l13y;P_Arm.l13z]).genProp;

EPF2EPCOM=System.Space.genEdge({'EPF2EPCOM' 'ExoPanFrame' 'ExoPanCOM'});
EL12F2EL1COM=System.Space.genEdge({'EL12F2EL1COM' 'ExoLink12Frame' 'ExoLink1COM'});
EL23F2EL2COM=System.Space.genEdge({'EL23F2EL2COM' 'ExoLink23Frame' 'ExoLink2COM'});
EL34F2EL3COM=System.Space.genEdge({'EL34F2EL3COM' 'ExoLink34Frame' 'ExoLink3COM'});
EL45F2EL4COM=System.Space.genEdge({'EL45F2EL4COM' 'ExoLink45Frame' 'ExoLink4COM'});
EAF2ELJCOM=System.Space.genEdge({'EAF2ELJCOM' 'ExoAttachFrame' 'ExoLinkJointCOM'});

EPF2EPCOM.setTransDis([P_COM2.lm22x;P_COM2.lm22y;P_COM2.lm22z]).genProp;
EL12F2EL1COM.setTransDis([P_COM2.lm23x;P_COM2.lm23y;P_COM2.lm23z]).genProp;
EL23F2EL2COM.setTransDis([P_COM2.lm24x;P_COM2.lm24y;P_COM2.lm24z]).genProp;
EL34F2EL3COM.setTransDis([P_COM2.lm25x;P_COM2.lm25y;P_COM2.lm25z]).genProp;
EL45F2EL4COM.setTransDis([P_COM2.lm26x;P_COM2.lm26y;P_COM2.lm26z]).genProp;
EAF2ELJCOM.setTransDis([P_COM2.lm27x;P_COM2.lm27y;P_COM2.lm27z]).genProp;
                 
%% Generate Kinematic System
System.Space.plotGraph(1);
System.Space.makeNumKinematics();