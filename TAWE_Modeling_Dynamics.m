qqdot=System.getContVector(1);


%% Declare Forearm Body
[ArmHandLoadBody,ArmHandBody]...
=System.genBody({
                'ArmHandLoadBody' 'ArmIMU2AdjustFrame';
                'ArmHandBody' 'ArmHandCOM';
                });

IArmHandLoadBody=diag([P_Uncertain.i12xx;P_Uncertain.i12yy;P_Uncertain.i12zz])...
                +[0,P_Uncertain.i12xy,P_Uncertain.i12xz;zeros(1,2),P_Uncertain.i12yz;zeros(1,3)]...
                +[0,P_Uncertain.i12xy,P_Uncertain.i12xz;zeros(1,2),P_Uncertain.i12yz;zeros(1,3)].';
ArmHandLoadBody.setProp(P_Uncertain.m12,IArmHandLoadBody);
ArmHandLoadBody.genForce({'forceArmHandLoadBodypt1' Handmpt1Frame WORLD}).setProp([0;0;-U_Uncertain.handmpt1*P_Sys.g]);
ArmHandLoadBody.genForce({'forceArmHandLoadBodypt2' Handmpt2Frame WORLD}).setProp([0;0;-U_Uncertain.handmpt2*P_Sys.g]);
ArmHandLoadBody.genForce({'forceArmHandLoadBodypt3' Handmpt3Frame WORLD}).setProp([0;0;-U_Uncertain.handmpt3*P_Sys.g]);
ArmHandLoadBody.genForce({'forceArmHandLoadBodypt4' Handmpt4Frame WORLD}).setProp([0;0;-U_Uncertain.handmpt4*P_Sys.g]);

IArmHandBody=diag([P_Inert1.i13xx;P_Inert1.i13yy;P_Inert1.i13zz])...
                +[0,P_Inert1.i13xy,P_Inert1.i13xz;zeros(1,2),P_Inert1.i13yz;zeros(1,3)]...
                +[0,P_Inert1.i13xy,P_Inert1.i13xz;zeros(1,2),P_Inert1.i13yz;zeros(1,3)].';
ArmHandBody.setProp(P_Inert1.m13,IArmHandBody).genForce({'GArmHandBody' ArmHandCOM WORLD}).setProp([0;0;-P_Inert1.m13*P_Sys.gUncertain]);

%% Declare Forearm Inputs
ArmHandBody.genTorque({'ArmDevTorque' ArmHandCOM ArmBaseDevFrame}).setProp([0;0;U1.u1dev]);
ArmHandBody.genTorque({'ArmFlexTorque' ArmHandCOM ArmDevFlexFrame}).setProp([U1.u1flex;0;0]);


%% Declare Forearm Damper
System.genDamper({'ArmDamper' 3 1}).setProp(P_Sys.b1act,qqdot(1:3));


%% Declare Exoskeleton Body
[ExoPanBody,ExoLink1Body,ExoLink2Body,ExoLink3Body,ExoLink4Body,ExoLinkJointBody]...
=System.genBody({
            'ExoPanBody' 'ExoPanCOM';
            'ExoLink1Body' 'ExoLink1COM';
            'ExoLink2Body' 'ExoLink2COM';
            'ExoLink3Body' 'ExoLink3COM';
            'ExoLink4Body' 'ExoLink4COM';
            'ExoLinkJointBody' 'ExoLinkJointCOM';
            });
       
IExoPanBody=diag([P_Inert2.i22xx;P_Inert2.i22yy;P_Inert2.i22zz])...
                +[0,P_Inert2.i22xy,P_Inert2.i22xz;zeros(1,2),P_Inert2.i22yz;zeros(1,3)]...
                +[0,P_Inert2.i22xy,P_Inert2.i22xz;zeros(1,2),P_Inert2.i22yz;zeros(1,3)].';
IExoLink1Body=diag([P_Inert2.i23xx;P_Inert2.i23yy;P_Inert2.i23zz])...
                +[0,P_Inert2.i23xy,P_Inert2.i23xz;zeros(1,2),P_Inert2.i23yz;zeros(1,3)]...
                +[0,P_Inert2.i23xy,P_Inert2.i23xz;zeros(1,2),P_Inert2.i23yz;zeros(1,3)].';
IExoLink2Body=diag([P_Inert2.i24xx;P_Inert2.i24yy;P_Inert2.i24zz])...
                +[0,P_Inert2.i24xy,P_Inert2.i24xz;zeros(1,2),P_Inert2.i24yz;zeros(1,3)]...
                +[0,P_Inert2.i24xy,P_Inert2.i24xz;zeros(1,2),P_Inert2.i24yz;zeros(1,3)].';
IExoLink3Body=diag([P_Inert2.i25xx;P_Inert2.i25yy;P_Inert2.i25zz])...
                +[0,P_Inert2.i25xy,P_Inert2.i25xz;zeros(1,2),P_Inert2.i25yz;zeros(1,3)]...
                +[0,P_Inert2.i25xy,P_Inert2.i25xz;zeros(1,2),P_Inert2.i25yz;zeros(1,3)].';
IExoLink4Body=diag([P_Inert2.i26xx;P_Inert2.i26yy;P_Inert2.i26zz])...
                +[0,P_Inert2.i26xy,P_Inert2.i26xz;zeros(1,2),P_Inert2.i26yz;zeros(1,3)]...
                +[0,P_Inert2.i26xy,P_Inert2.i26xz;zeros(1,2),P_Inert2.i26yz;zeros(1,3)].';
IExoLinkJointBody=diag([P_Inert2.i27xx;P_Inert2.i27yy;P_Inert2.i27zz])...
                +[0,P_Inert2.i27xy,P_Inert2.i27xz;zeros(1,2),P_Inert2.i27yz;zeros(1,3)]...
                +[0,P_Inert2.i27xy,P_Inert2.i27xz;zeros(1,2),P_Inert2.i27yz;zeros(1,3)].';
            
ExoPanBody.setProp(P_Inert2.m22,IExoPanBody).genForce({'GExoPanBody' ExoPanCOM WORLD}).setProp([0;0;-P_Inert2.m22*P_Sys.g]);
ExoLink1Body.setProp(P_Inert2.m23,IExoLink1Body).genForce({'GExoLink1Body' ExoLink1COM WORLD}).setProp([0;0;-P_Inert2.m23*P_Sys.g]);
ExoLink2Body.setProp(P_Inert2.m24,IExoLink2Body).genForce({'GExoLink2Body' ExoLink2COM WORLD}).setProp([0;0;-P_Inert2.m24*P_Sys.g]);
ExoLink3Body.setProp(P_Inert2.m25,IExoLink3Body).genForce({'GExoLink3Body' ExoLink3COM WORLD}).setProp([0;0;-P_Inert2.m25*P_Sys.g]);
ExoLink4Body.setProp(P_Inert2.m26,IExoLink4Body).genForce({'GExoLink4Body' ExoLink4COM WORLD}).setProp([0;0;-P_Inert2.m26*P_Sys.g]);
ExoLinkJointBody.setProp(P_Inert2.m27,IExoLinkJointBody).genForce({'GExoLinkJointBody' ExoLinkJointCOM WORLD}).setProp([0;0;-P_Inert2.m27*P_Sys.g]);


%% Declare Exoskeleton Inputs
ExoPanBody.genTorque({'ExoPanTorque' ExoPanCOM ExoPanFrame}).setProp([0;0;U2.u21]);
ExoPanBody.genTorque({'ExoCounterLink1Torque' ExoPanCOM ExoLink12Frame}).setProp([-U2.u22;0;0]);
ExoLink1Body.genTorque({'ExoLink1Torque' ExoLink1COM ExoLink12Frame}).setProp([U2.u22;0;0]);


%% Declare Exoskeleton Damper
System.genDamper({'ExoActDamper' 6 1}).setProp(P_Sys.b2act,qqdot(4:end));


%% Declare Constraint Attach Mode
[AttachExoEndDx,AttachExoEndDy,AttachExoEndDz,...
AttachExoEndRx,AttachExoEndRy,AttachExoEndRz,SupConsY]...
=System.genCons({
                 'AttachExoEndDx' 'l_AEEdx';
                 'AttachExoEndDy' 'l_AEEdy';
                 'AttachExoEndDz' 'l_AEEdz';
                 'AttachExoEndRx' 'l_AEErx';
                 'AttachExoEndRy' 'l_AEEry';
                 'AttachExoEndRz' 'l_AEErz';
                 'SupConsY' 'l_scy';
               });
ArmAttachEnd=[eye(3),zeros(3,1)]*(ArmIMU2Frame.fwdTF(ArmIMU1Frame))*[zeros(3,1);1];
ExoAttachEnd=[eye(3),zeros(3,1)]*(ExoIMU2Frame.fwdTF(ExoIMU1Frame))*[zeros(3,1);1];
AEEPointExpr=ArmAttachEnd-ExoAttachEnd;
AttachExoEndDx.setHCons(AEEPointExpr(1),P_Sys.c);
AttachExoEndDy.setHCons(AEEPointExpr(2),P_Sys.c);
AttachExoEndDz.setHCons(AEEPointExpr(3),P_Sys.c);
AttachArmRotQuat=quatMultiply(quatConj([1;0;0;0]),ArmIMU2Frame.getAngSym(0));
AttachExoRotQuat=quatMultiply(quatConj([1;0;0;0]),ExoIMU2Frame.getAngSym(0));
AEEAttitudeExpr=[zeros(3,1),eye(3)]*quatMultiply(AttachArmRotQuat,quatConj(AttachExoRotQuat));
AttachExoEndRx.setHCons(AEEAttitudeExpr(1),P_Sys.c);
AttachExoEndRy.setHCons(AEEAttitudeExpr(2),P_Sys.c);
AttachExoEndRz.setHCons(AEEAttitudeExpr(3),P_Sys.c);

rotMat=ypr2Mat(kappa);
rotQuat=ypr2Quat(kappa);
halfCos=sqrt((rotQuat(1)+1)/2);
rotQuatSquareRoot=[halfCos;rotQuat(2:4)/2/halfCos];
vel3=pDiff(quat2YPR(rotQuatSquareRoot)*2,t);
xx3.setExpr(jSimplify(vel3(1)));
yy3.setExpr(jSimplify(vel3(2)));
zz3.setExpr(jSimplify(vel3(3)));
lockyExpr=[1 0 0 0]*quatMultiply([0;[0;1;0]],rotQuat) + sin(xx3.dot(0)/2)*sin(zz3.dot(0)/2)/2;
SupConsY.setHCons(lockyExpr,P_Sys.c);

