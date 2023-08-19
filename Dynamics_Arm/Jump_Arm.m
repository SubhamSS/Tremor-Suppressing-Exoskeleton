function [flow,flowTraj,q,p,s]=Jump_Arm(flow,t,q,p,u,s,l)
%% Author: Jiamin Wang; Updated: 2021-12-15;
    act=0;
    count=0;
    flowTraj=zeros(100+1,1);
    while(1)
        flowTraj(count+1)=flow;
        act=0;
%SWITCHCASE_

	switch flow
	    case 1
	        %No Jump Here
	end

        if(act==0)
            break;
        end
        if(count>100)
            flow=-flow;%Enter zeno or other infinite jump problem, quit and give an error.
            break;
        end
    end

	function [TF,Vel,Cor,Jac,TransDis,Quat,MM,GFI,GFF,GFU,JacU,JacCons,CorCons,GFCons,CenMat,CorMatLeft,CorMatRight]=System_Arm(t,q,p,u,s)
	%% Author: Jiamin Wang; Updated: 2021-12-15;
	    
	    [TF,Vel,Cor,Jac,TransDis,Quat]=numKinematics_Arm(t,q,p,u,s);
	    [MM,GFI,CenMat,CorMatLeft,CorMatRight]=Inertial_Arm(t,q,p,u,s,TF,Vel,Cor,Jac,true);
	    [GFF]=Force_Arm(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
	    [GFU,JacU]=Input_Arm(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
	    [JacCons,CorCons,GFCons]=Constraint_Arm(t,q,p,u,s,TransDis,Vel,Cor,Jac,Quat);
		function [frameTF,frameVel,frameCorAcc,frameJacobian,frameTransDis,frameRotQuat]=numKinematics_Arm(t,q,p,u,s)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    LinkTF=linkTF(t,q,p,u,s);
		    LinkVel=linkVel(t,q,p,u,s);
		    LinkCorAcc=linkCorAcc(t,q,p,u,s);
		    LinkJacobian=linkJacobian(t,q,p,u,s);
		    LinkRotQuat=linkRotQuat(t,q,p,u,s);
		    FrameNum=11;
		    DOF=numel(q)/2;
		    
		    LinkTF=reshape((LinkTF),[4,4,FrameNum]);
		    LinkJacobian=reshape((LinkJacobian),[6,DOF,FrameNum]);
		    
		    frameTF=zeros(4,4,FrameNum);
		    frameTransDis=zeros(3,FrameNum);
		    frameVel=zeros(6,FrameNum);
		    frameCorAcc=zeros(6,FrameNum);
		    frameJacobian=zeros(6,DOF,FrameNum);
		    frameRotQuat=zeros(4,FrameNum);
		    
		    frameTF(:,:,1)=LinkTF(:,:,1);
		    frameTransDis(:,1)=frameTF(1:3,4,1);
		    frameRotQuat(:,1)=LinkRotQuat(:,1);
		    frameVel(:,1)=LinkVel(:,1);
		    frameCorAcc(:,1)=LinkCorAcc(:,1);
		    frameJacobian(:,:,1)=LinkJacobian(:,:,1);
		    frameCheckList=zeros(1,FrameNum);
		    frameCheckList(1)=1;
		    for frame=1:FrameNum
		        framePath=0;
		%SWITCHCASE_
				switch frame
				    case 1
				        framePath=[1];
				    case 2
				        framePath=[1  2];
				    case 3
				        framePath=[1  2  3];
				    case 4
				        framePath=[1  2  3  4];
				    case 5
				        framePath=[1  2  3  4  5];
				    case 6
				        framePath=[1  2  3  4  5  6];
				    case 7
				        framePath=[1  2  3  4  5  6  7];
				    case 8
				        framePath=[1  2  3  4  5  6  8];
				    case 9
				        framePath=[1  2  3  4  5  6  9];
				    case 10
				        framePath=[1   2   3   4   5   6  10];
				    case 11
				        framePath=[1   2   3   4   5   6  11];
				end
		        [frameCheckList,LinkTF,frameTF,LinkVel,frameVel,LinkCorAcc,frameCorAcc,LinkJacobian,frameJacobian,frameTransDis,LinkRotQuat,frameRotQuat]=frameKineRecur(framePath,frameCheckList,LinkTF,frameTF,LinkVel,frameVel,LinkCorAcc,frameCorAcc,LinkJacobian,frameJacobian,frameTransDis,LinkRotQuat,frameRotQuat);
		    end
		    function output_=localCrossFunc(aa_,bb_)
		        output_=[aa_(2) * bb_(3) - aa_(3) * bb_(2); aa_(3) * bb_(1) - aa_(1) * bb_(3) ; aa_(1) * bb_(2) - aa_(2) * bb_(1) ];
		    end
		    function output_=skew3(w_)
		        output_=[0 -w_(3) w_(2) ; w_(3) 0 -w_(1) ; -w_(2) w_(1) 0 ];
		    end
		    function p_ = quatMultiply(q_, r_)
		        p_ = [q_(1).*r_(1) - q_(2).*r_(2) - q_(3).*r_(3) - q_(4).*r_(4);
		             q_(1).*r_(2) + r_(1).*q_(2) + q_(3).*r_(4) - q_(4).*r_(3);
		             q_(1).*r_(3) + r_(1).*q_(3) + q_(4).*r_(2) - q_(2).*r_(4);
		             q_(1).*r_(4) + r_(1).*q_(4) + q_(2).*r_(3) - q_(3).*r_(2)];
		    end
		    function [frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_]=frameKineRecur(framePath_,frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_)
		        curNum_=framePath_(end);
		        if(~frameCheckList_(curNum_))
		            preNum_=framePath_(end-1);
		            if(frameCheckList_(preNum_))
		                frameTF_(:,:,curNum_)=frameTF_(:,:,preNum_)*TF_(:,:,curNum_);
		                frameTransDis_(:,curNum_)=frameTF_(1:3,4,curNum_);
		                frameRotQuat_(:,curNum_)=quatMultiply(frameRotQuat_(:,preNum_),RotQuat_(:,curNum_));
		                frameRotQuat_(:,curNum_)=frameRotQuat_(:,curNum_)./sqrt(frameRotQuat_(1,curNum_)^2 + ...
		                                         frameRotQuat_(2,curNum_)^2 + frameRotQuat_(3,curNum_)^2 + frameRotQuat_(4,curNum_)^2);
		                curTransVel=Vel_(1:3,curNum_);
		                curAngVel=Vel_(4:6,curNum_);
		                curCorTransAcc=CorAcc_(1:3,curNum_);
		                curCorAngAcc=CorAcc_(4:6,curNum_);
		                curTransJacobian=Jacobian_(1:3,:,curNum_);
		                curAngJacobian=Jacobian_(4:6,:,curNum_);
		                preTransVel=frameVel_(1:3,preNum_);
		                preAngVel=frameVel_(4:6,preNum_);
		                preCorTransAcc=frameCorAcc_(1:3,preNum_);
		                preCorAngAcc=frameCorAcc_(4:6,preNum_);
		                preTransJacobian=frameJacobian_(1:3,:,preNum_);
		                preAngJacobian=frameJacobian_(4:6,:,preNum_);
		                
		                preRotMat=frameTF_(1:3,1:3,preNum_);
		                preRotCurTransDis=preRotMat*TF_(1:3,4,curNum_);
		                
		                frameVel_(1:3,curNum_)=preTransVel+preRotMat*curTransVel+localCrossFunc(preAngVel,preRotCurTransDis);
		                frameVel_(4:6,curNum_)=preAngVel+preRotMat*curAngVel;
		                frameCorAcc_(1:3,curNum_)=preCorTransAcc+localCrossFunc(preAngVel,localCrossFunc(preAngVel,preRotCurTransDis))...
		                                        +preRotMat*curCorTransAcc+localCrossFunc(preCorAngAcc,preRotCurTransDis)...
		                                        +2*localCrossFunc(preAngVel,preRotMat*curTransVel);
		                frameCorAcc_(4:6,curNum_)=preCorAngAcc+preRotMat*curCorAngAcc+localCrossFunc(preAngVel,preRotMat*curAngVel);
		                %skewTransDis=[0 -preRotCurTransDis(3) preRotCurTransDis(2) ; preRotCurTransDis(3) 0 -preRotCurTransDis(1) ; -preRotCurTransDis(2) preRotCurTransDis(1) 0 ];
		                skewTransDis=skew3(preRotCurTransDis);
		                frameJacobian_(1:3,:,curNum_)=preTransJacobian+preRotMat*curTransJacobian-skewTransDis*preAngJacobian;
		                frameJacobian_(4:6,:,curNum_)=preAngJacobian+preRotMat*curAngJacobian;
		                frameCheckList_(curNum_)=1;
		            else
		                [frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_]=frameKineRecur(framePath_(1:end-1),frameCheckList_,TF_,frameTF_,Vel_,frameVel_,CorAcc_,frameCorAcc_,Jacobian_,frameJacobian_,frameTransDis_,RotQuat_,frameRotQuat_);
		            end
		        end
		    end
			function out1 = linkTF(t,in2,in3,in4,in5)
			%linkTF
			%    OUT1 = linkTF(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:11:31
			l11x = in3(7,:);
			l11y = in3(8,:);
			l12x = in3(10,:);
			l11z = in3(9,:);
			l12y = in3(11,:);
			l12z = in3(12,:);
			lm13x = in3(25,:);
			lm13y = in3(26,:);
			lm13z = in3(27,:);
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			rTetra = in3(6,:);
			rf1IMU1x = in3(16,:);
			rf1IMU1y = in3(17,:);
			rf1IMU2x = in3(22,:);
			rf1IMU1z = in3(18,:);
			rf1IMU2y = in3(23,:);
			rf1IMU2z = in3(24,:);
			rf1Devx = in3(19,:);
			rf1Devy = in3(20,:);
			rf1Devz = in3(21,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf1IMU1x);
			t6 = cos(rf1IMU1y);
			t7 = cos(rf1IMU2x);
			t8 = cos(rf1IMU1z);
			t9 = cos(rf1IMU2y);
			t10 = cos(rf1IMU2z);
			t11 = cos(rf1Devx);
			t12 = cos(rf1Devy);
			t13 = cos(rf1Devz);
			t14 = sin(q1dev__dt_0_);
			t15 = sin(q1flex__dt_0_);
			t16 = sin(q1sup__dt_0_);
			t17 = sin(rf1IMU1x);
			t18 = sin(rf1IMU1y);
			t19 = sin(rf1IMU2x);
			t20 = sin(rf1IMU1z);
			t21 = sin(rf1IMU2y);
			t22 = sin(rf1IMU2z);
			t23 = sin(rf1Devx);
			t24 = sin(rf1Devy);
			t25 = sin(rf1Devz);
			t26 = sqrt(3.0);
			t27 = sqrt(6.0);
			t28 = q1dev__dt_0_./2.0;
			t29 = q1flex__dt_0_./2.0;
			t30 = q1sup__dt_0_./2.0;
			t31 = rf1IMU1x./2.0;
			t32 = rf1IMU1y./2.0;
			t33 = rf1IMU2x./2.0;
			t34 = rf1IMU1z./2.0;
			t35 = rf1IMU2y./2.0;
			t36 = rf1IMU2z./2.0;
			t37 = rf1Devx./2.0;
			t38 = rf1Devy./2.0;
			t39 = rf1Devz./2.0;
			t40 = cos(t28);
			t41 = cos(t29);
			t42 = cos(t31);
			t43 = cos(t33);
			t44 = cos(t34);
			t45 = cos(t36);
			t46 = cos(t37);
			t47 = cos(t39);
			t56 = (rTetra.*t26)./3.0;
			t57 = (rTetra.*t27)./6.0;
			t48 = t40.^2;
			t49 = t41.^2;
			t50 = t42.^2;
			t51 = t43.^2;
			t52 = t44.^2;
			t53 = t45.^2;
			t54 = t46.^2;
			t55 = t47.^2;
			t58 = -t56;
			t59 = -t57;
			mt1 = [1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,t6.*t8,t6.*t20,-t18,0.0,-t5.*t20+t8.*t17.*t18,t50.*-2.0-t52.*2.0+t50.*t52.*4.0+t42.*t44.*cos(t32).*sin(t31).*sin(t32).*sin(t34).*8.0+1.0,t6.*t17,0.0,t17.*t20+t5.*t8.*t18,-t8.*t17+t5.*t18.*t20,t5.*t6,0.0,0.0,0.0,0.0,1.0,t12.*t13,t12.*t25,-t24,0.0,-t11.*t25+t13.*t23.*t24,t54.*-2.0-t55.*2.0+t54.*t55.*4.0+t46.*t47.*cos(t38).*sin(t37).*sin(t38).*sin(t39).*8.0+1.0,t12.*t23,0.0,t23.*t25+t11.*t13.*t24,-t13.*t23+t11.*t24.*t25,t11.*t12,0.0,l11x,l11y,l11z,1.0,t2.*t4,t4.*t14,-t16,0.0,-t3.*t14+t2.*t15.*t16];
			mt2 = [t48.*-2.0-t49.*2.0+t48.*t49.*4.0+t40.*t41.*cos(t30).*sin(t28).*sin(t29).*sin(t30).*8.0+1.0,t4.*t15,0.0,t14.*t15+t2.*t3.*t16,-t2.*t15+t3.*t14.*t16,t3.*t4,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,l12x,l12y,l12z,1.0,t9.*t10,t9.*t22,-t21,0.0,-t7.*t22+t10.*t19.*t21,t51.*-2.0-t53.*2.0+t51.*t53.*4.0+t43.*t45.*cos(t35).*sin(t33).*sin(t35).*sin(t36).*8.0+1.0,t9.*t19,0.0,t19.*t22+t7.*t10.*t21,-t10.*t19+t7.*t21.*t22,t7.*t9,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm13x,lm13y,lm13z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,-rTetra,t59,t58,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,rTetra,t59];
			mt3 = [t58,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,t59,rTetra.*t26.*(2.0./3.0),1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,(rTetra.*t27)./2.0,0.0,1.0];
			out1 = reshape([mt1,mt2,mt3],4,44);
			end
			function out1 = linkVel(t,in2,in3,in4,in5)
			%linkVel
			%    OUT1 = linkVel(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:11:30
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(4,:);
			q1flex__dt_1_ = in2(5,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(6,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = sin(q1dev__dt_0_);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-q1sup__dt_1_.*t4+q1flex__dt_1_.*t2.*t3,q1sup__dt_1_.*t2+q1flex__dt_1_.*t3.*t4,q1dev__dt_1_-q1flex__dt_1_.*sin(q1sup__dt_0_),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,11]);
			end
			function out1 = linkRotQuat(t,in2,in3,in4,in5)
			%linkRotQuat
			%    OUT1 = linkRotQuat(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:11:30
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			rf1IMU1x = in3(16,:);
			rf1IMU1y = in3(17,:);
			rf1IMU2x = in3(22,:);
			rf1IMU1z = in3(18,:);
			rf1IMU2y = in3(23,:);
			rf1IMU2z = in3(24,:);
			rf1Devx = in3(19,:);
			rf1Devy = in3(20,:);
			rf1Devz = in3(21,:);
			t2 = q1dev__dt_0_./2.0;
			t3 = q1flex__dt_0_./2.0;
			t4 = q1sup__dt_0_./2.0;
			t5 = rf1IMU1x./2.0;
			t6 = rf1IMU1y./2.0;
			t7 = rf1IMU2x./2.0;
			t8 = rf1IMU1z./2.0;
			t9 = rf1IMU2y./2.0;
			t10 = rf1IMU2z./2.0;
			t11 = rf1Devx./2.0;
			t12 = rf1Devy./2.0;
			t13 = rf1Devz./2.0;
			t14 = cos(t2);
			t15 = cos(t3);
			t16 = cos(t4);
			t17 = cos(t5);
			t18 = cos(t6);
			t19 = cos(t7);
			t20 = cos(t8);
			t21 = cos(t9);
			t22 = cos(t10);
			t23 = cos(t11);
			t24 = cos(t12);
			t25 = cos(t13);
			t26 = sin(t2);
			t27 = sin(t3);
			t28 = sin(t4);
			t29 = sin(t5);
			t30 = sin(t6);
			t31 = sin(t7);
			t32 = sin(t8);
			t33 = sin(t9);
			t34 = sin(t10);
			t35 = sin(t11);
			t36 = sin(t12);
			t37 = sin(t13);
			out1 = reshape([1.0,0.0,0.0,0.0,t17.*t18.*t20+t29.*t30.*t32,t18.*t20.*t29-t17.*t30.*t32,t17.*t20.*t30+t18.*t29.*t32,t17.*t18.*t32-t20.*t29.*t30,t23.*t24.*t25+t35.*t36.*t37,t24.*t25.*t35-t23.*t36.*t37,t23.*t25.*t36+t24.*t35.*t37,t23.*t24.*t37-t25.*t35.*t36,t14.*t15.*t16+t26.*t27.*t28,t14.*t16.*t27-t15.*t26.*t28,t14.*t15.*t28+t16.*t26.*t27,t15.*t16.*t26-t14.*t27.*t28,1.0,0.0,0.0,0.0,t19.*t21.*t22+t31.*t33.*t34,t21.*t22.*t31-t19.*t33.*t34,t19.*t22.*t33+t21.*t31.*t34,t19.*t21.*t34-t22.*t31.*t33,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[4,11]);
			end
			function out1 = linkCorAcc(t,in2,in3,in4,in5)
			%linkCorAcc
			%    OUT1 = linkCorAcc(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:11:30
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(4,:);
			q1flex__dt_1_ = in2(5,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(6,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = sin(q1dev__dt_0_);
			t5 = sin(q1sup__dt_0_);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-q1dev__dt_1_.*q1sup__dt_1_.*t2-q1dev__dt_1_.*q1flex__dt_1_.*t3.*t4-q1flex__dt_1_.*q1sup__dt_1_.*t2.*t5,-q1dev__dt_1_.*q1sup__dt_1_.*t4+q1dev__dt_1_.*q1flex__dt_1_.*t2.*t3-q1flex__dt_1_.*q1sup__dt_1_.*t4.*t5,-q1flex__dt_1_.*q1sup__dt_1_.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,11]);
			end
			function out1 = linkJacobian(t,in2,in3,in4,in5)
			%linkJacobian
			%    OUT1 = linkJacobian(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:11:31
			q1dev__dt_0_ = in2(1,:);
			q1sup__dt_0_ = in2(3,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = sin(q1dev__dt_0_);
			mt1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t2.*t3,t3.*t4,-sin(q1sup__dt_0_),0.0,0.0,0.0,-t4,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			out1 = reshape([mt1,mt2],6,33);
			end
		end
		function [inertM,inertGF,inertCenMat,inertCorMatLeft,inertCorMatRight]=Inertial_Arm(t,q,p,u,s,TF_Global,Vel_Global,Cor_Global,Jac_Global,Full_Flag_)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    BaseFrameList=[6  7];
		    MassList=Mass_Arm(t,q,p,s,u);
		    if Full_Flag_
		        MomentList=reshape(Moment_Arm(t,q,p,s,u),[3 3 numel(MassList)]);
		    end
		    
		    inertCorMatLeft=zeros(numel(q)/2,6*numel(MassList));
		    inertCorMatRight=zeros(6*numel(MassList),numel(q)/2);
		    inertM=zeros(numel(q)/2,numel(q)/2,numel(MassList));
		    inertGF=zeros(numel(q)/2,numel(MassList));
		    inertCenMat=zeros(numel(q)/2,numel(q)/2,numel(MassList));
		    for bodyNum=1:numel(MassList)
		        if Full_Flag_
		            rotor=TF_Global(1:3,1:3,BaseFrameList(bodyNum));
		            m=MassList(bodyNum);
		            I=rotor*MomentList(:,:,bodyNum)*rotor.';
		            w=Vel_Global(4:6,BaseFrameList(bodyNum));
		            a=Cor_Global(1:3,BaseFrameList(bodyNum));
		            alpha=Cor_Global(4:6,BaseFrameList(bodyNum));
		        end
		        vJacobian=Jac_Global(1:3,:,BaseFrameList(bodyNum));
		        wJacobian=Jac_Global(4:6,:,BaseFrameList(bodyNum));
		        
		        inertCorMatRight(6*bodyNum+(-5:0),:)=[vJacobian;wJacobian];
		        if Full_Flag_
		            % thisInertCorMatLeft=[m*vJacobian;I*wJacobian].';
		            % thisInertCorMatRight=[vJacobian;wJacobian]; %Take the time derivative of this to calculate coriolis numerically
		            inertCorMatLeft(:,6*bodyNum+(-5:0))=[m*vJacobian;I*wJacobian].';
		            inertCenMat(:,:,bodyNum)=- wJacobian.'* skew3_InertDynamicsOnly_(I*w) * wJacobian; % On the same side of Mqddot
		            % inertM(:,:,bodyNum)=thisInertCorMatLeft*thisInertCorMatRight;
		            % inertGF(:,bodyNum)= - thisInertCorMatLeft*[a;alpha] - wJacobian.'* cross(w,I*w); % on the different side of Mqddot
		            inertM(:,:,bodyNum)=vJacobian.'*m*vJacobian+wJacobian.'*I*wJacobian;
		            inertGF(:,bodyNum)=- vJacobian.'*m*a - wJacobian.'*(I*alpha+cross(w,I*w)); % on the different side of Mqddot
		        end
		    end
		    function output_=skew3_InertDynamicsOnly_(w_)
		        output_=[0 -w_(3) w_(2) ; w_(3) 0 -w_(1) ; -w_(2) w_(1) 0 ];
		    end
			function out1 = Mass_Arm(t,in2,in3,in4,in5)
			%Mass_Arm
			%    OUT1 = Mass_Arm(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:05
			m12 = in3(35,:);
			m13 = in3(28,:);
			out1 = [m12,m13];
			end
			function out1 = Moment_Arm(t,in2,in3,in4,in5)
			%Moment_Arm
			%    OUT1 = Moment_Arm(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:05
			i12xx = in3(36,:);
			i12xy = in3(37,:);
			i13xx = in3(29,:);
			i12xz = in3(38,:);
			i12yy = in3(39,:);
			i13xy = in3(30,:);
			i12yz = in3(40,:);
			i13xz = in3(31,:);
			i13yy = in3(32,:);
			i12zz = in3(41,:);
			i13yz = in3(33,:);
			i13zz = in3(34,:);
			out1 = reshape([i12xx,i12xy,i12xz,i12xy,i12yy,i12yz,i12xz,i12yz,i12zz,i13xx,i13xy,i13xz,i13xy,i13yy,i13yz,i13xz,i13yz,i13zz],[3,6]);
			end
		end
		function [CollectGF]=Force_Arm(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		 
		    CollectGF=zeros(numel(q)/2,1);
		    ForceNum=2;
		    if ForceNum>0
		        CollectGF=zeros(numel(q)/2,ForceNum);
		    end
		    SubGF=zeros(numel(q)/2,1);
		    for FCount=1:ForceNum
		%SWITCHCASE_
				switch FCount
				    case 1
				        SubFrame=[0];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ForceFrameVecJac=ForceFVJac_1(t,q,p,u,s,SubSubs);
				        SubJac=ForceJac_1(t,q,p,u,s,SubSubs);
				        SubEff=ForceEff_1(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubJac=SubJac+ForceFrameVecJac(:,:,sfCount)*JacSubs(:,:,sfCount);
				            end
				        end
				        SubGF=SubJac.'*SubEff;
				    case 2
				        SubFrame=[7];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[7];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_2(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				end
		        CollectGF(:,FCount)=SubGF;
		    end
			function out1 = ForceJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceJac_1
			%    OUT1 = ForceJac_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:05
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,3]);
			end
			function out1 = ForceJac_2(t,in2,in3,in4,in5,in6)
			%ForceJac_2
			%    OUT1 = ForceJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			out1 = 0.0;
			end
			function out1 = ForceEff_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceEff_1
			%    OUT1 = ForceEff_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			b1act = in3(4,:);
			q1dev__dt_1_ = in2(4,:);
			q1flex__dt_1_ = in2(5,:);
			q1sup__dt_1_ = in2(6,:);
			out1 = [-b1act.*q1dev__dt_1_;-b1act.*q1flex__dt_1_;-b1act.*q1sup__dt_1_];
			end
			function out1 = ForceEff_2(t,in2,in3,in4,in5,in6)
			%ForceEff_2
			%    OUT1 = ForceEff_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			gUncertain = in3(2,:);
			m13 = in3(28,:);
			out1 = [0.0;0.0;-gUncertain.*m13];
			end
			function out1 = ForceFVJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceFVJac_1
			%    OUT1 = ForceFVJac_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			out1 = 0.0;
			end
			function out1 = ForceFVJac_2(t,in2,in3,in4,in5,in6)
			%ForceFVJac_2
			%    OUT1 = ForceFVJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
		end
		function [CollectGF,CollectInputJac]=Input_Arm(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    CollectGF=zeros(numel(q)/2,1);
		    CollectInputJac=zeros(numel(q)/2,numel(u));
		    ForceNum=6;
		    if ForceNum>0
		        CollectGF=zeros(numel(q)/2,ForceNum);
		        CollectInputJac=zeros(numel(q)/2,numel(u),ForceNum);
		    end
		    SubGF=zeros(numel(q)/2,1);
		    SubInputJac=zeros(numel(q)/2,numel(u));
		    for FCount=1:ForceNum
		%SWITCHCASE_
				switch FCount
				    case 1
				        SubFrame=[8];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[8];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_1(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_1(t,q,p,u,s,SubSubs);
				    case 2
				        SubFrame=[9];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[9];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_2(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_2(t,q,p,u,s,SubSubs);
				    case 3
				        SubFrame=[10];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[10];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_3(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_3(t,q,p,u,s,SubSubs);
				    case 4
				        SubFrame=[11];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[11];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_4(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_4(t,q,p,u,s,SubSubs);
				    case 5
				        SubFrame=[7];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[-7];
				        RefFrame=[3];
				        SubJac=Jac_Global(4:6,:,-ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_5(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_5(t,q,p,u,s,SubSubs);
				    case 6
				        SubFrame=[7];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[-7];
				        RefFrame=[4];
				        SubJac=Jac_Global(4:6,:,-ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_6(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_6(t,q,p,u,s,SubSubs);
				end
		        CollectGF(:,FCount)=SubGF;
		        CollectInputJac(:,:,FCount)=SubInputJac;
		    end
			function out1 = ForceJac_1(t,in2,in3,in4,in5,in6)
			%ForceJac_1
			%    OUT1 = ForceJac_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			out1 = 0.0;
			end
			function out1 = ForceJac_2(t,in2,in3,in4,in5,in6)
			%ForceJac_2
			%    OUT1 = ForceJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			out1 = 0.0;
			end
			function out1 = ForceJac_3(t,in2,in3,in4,in5,in6)
			%ForceJac_3
			%    OUT1 = ForceJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			out1 = 0.0;
			end
			function out1 = ForceJac_4(t,in2,in3,in4,in5,in6)
			%ForceJac_4
			%    OUT1 = ForceJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			out1 = 0.0;
			end
			function out1 = ForceJac_5(t,in2,in3,in4,in5,in6)
			%ForceJac_5
			%    OUT1 = ForceJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			out1 = 0.0;
			end
			function out1 = ForceJac_6(t,in2,in3,in4,in5,in6)
			%ForceJac_6
			%    OUT1 = ForceJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:09
			out1 = 0.0;
			end
			function out1 = InputEff_1(t,in2,in3,in4,in5,in6)
			%InputEff_1
			%    OUT1 = InputEff_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:06
			g = in3(1,:);
			handmpt1 = in4(3,:);
			out1 = [0.0;0.0;-g.*handmpt1];
			end
			function out1 = InputEff_2(t,in2,in3,in4,in5,in6)
			%InputEff_2
			%    OUT1 = InputEff_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			g = in3(1,:);
			handmpt2 = in4(4,:);
			out1 = [0.0;0.0;-g.*handmpt2];
			end
			function out1 = InputEff_3(t,in2,in3,in4,in5,in6)
			%InputEff_3
			%    OUT1 = InputEff_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			g = in3(1,:);
			handmpt3 = in4(5,:);
			out1 = [0.0;0.0;-g.*handmpt3];
			end
			function out1 = InputEff_4(t,in2,in3,in4,in5,in6)
			%InputEff_4
			%    OUT1 = InputEff_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			g = in3(1,:);
			handmpt4 = in4(6,:);
			out1 = [0.0;0.0;-g.*handmpt4];
			end
			function out1 = InputEff_5(t,in2,in3,in4,in5,in6)
			%InputEff_5
			%    OUT1 = InputEff_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			u1dev = in4(1,:);
			out1 = [0.0;0.0;u1dev];
			end
			function out1 = InputEff_6(t,in2,in3,in4,in5,in6)
			%InputEff_6
			%    OUT1 = InputEff_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:09
			u1flex = in4(2,:);
			out1 = [u1flex;0.0;0.0];
			end
			function out1 = InputEffJac_1(t,in2,in3,in4,in5,in6)
			%InputEffJac_1
			%    OUT1 = InputEffJac_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputEffJac_2(t,in2,in3,in4,in5,in6)
			%InputEffJac_2
			%    OUT1 = InputEffJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputEffJac_3(t,in2,in3,in4,in5,in6)
			%InputEffJac_3
			%    OUT1 = InputEffJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputEffJac_4(t,in2,in3,in4,in5,in6)
			%InputEffJac_4
			%    OUT1 = InputEffJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g],[3,6]);
			end
			function out1 = InputEffJac_5(t,in2,in3,in4,in5,in6)
			%InputEffJac_5
			%    OUT1 = InputEffJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:09
			rf1IMU1x = in3(16,:);
			rf1IMU1y = in3(17,:);
			rf1IMU1z = in3(18,:);
			rf1Devx = in3(19,:);
			rf1Devy = in3(20,:);
			rf1Devz = in3(21,:);
			t2 = cos(rf1IMU1x);
			t3 = cos(rf1IMU1y);
			t4 = cos(rf1IMU1z);
			t5 = cos(rf1Devx);
			t6 = cos(rf1Devy);
			t7 = cos(rf1Devz);
			t8 = sin(rf1IMU1x);
			t9 = sin(rf1IMU1y);
			t10 = sin(rf1IMU1z);
			t11 = sin(rf1Devx);
			t12 = sin(rf1Devy);
			t13 = sin(rf1Devz);
			t16 = rf1IMU1x./2.0;
			t17 = rf1IMU1y./2.0;
			t18 = rf1IMU1z./2.0;
			t14 = t7.*t11;
			t15 = t11.*t13;
			t19 = cos(t16);
			t20 = cos(t18);
			t21 = t5.*t12.*t13;
			t22 = t5.*t7.*t12;
			t23 = t19.^2;
			t24 = t20.^2;
			t25 = -t21;
			t26 = t15+t22;
			t27 = t14+t25;
			out1 = reshape([t27.*(t2.*t10-t4.*t8.*t9)+t5.*t6.*(t8.*t10+t2.*t4.*t9)+t3.*t4.*t26,-t27.*(t23.*-2.0-t24.*2.0+t23.*t24.*4.0+t19.*t20.*cos(t17).*sin(t16).*sin(t17).*sin(t18).*8.0+1.0)-t5.*t6.*(t4.*t8-t2.*t9.*t10)+t3.*t10.*t26,-t9.*t26-t3.*t8.*t27+t2.*t3.*t5.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputEffJac_6(t,in2,in3,in4,in5,in6)
			%InputEffJac_6
			%    OUT1 = InputEffJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:09
			q1dev__dt_0_ = in2(1,:);
			q1sup__dt_0_ = in2(3,:);
			rf1IMU1x = in3(16,:);
			rf1IMU1y = in3(17,:);
			rf1IMU1z = in3(18,:);
			rf1Devx = in3(19,:);
			rf1Devy = in3(20,:);
			rf1Devz = in3(21,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = cos(rf1IMU1x);
			t5 = cos(rf1IMU1y);
			t6 = cos(rf1IMU1z);
			t7 = cos(rf1Devx);
			t8 = cos(rf1Devy);
			t9 = cos(rf1Devz);
			t10 = sin(q1dev__dt_0_);
			t11 = sin(q1sup__dt_0_);
			t12 = sin(rf1IMU1x);
			t13 = sin(rf1IMU1y);
			t14 = sin(rf1IMU1z);
			t15 = sin(rf1Devx);
			t16 = sin(rf1Devy);
			t17 = sin(rf1Devz);
			t21 = rf1IMU1x./2.0;
			t22 = rf1IMU1y./2.0;
			t23 = rf1IMU1z./2.0;
			t24 = rf1Devx./2.0;
			t25 = rf1Devy./2.0;
			t26 = rf1Devz./2.0;
			t18 = t7.*t17;
			t19 = t9.*t15;
			t20 = t15.*t17;
			t27 = cos(t21);
			t28 = cos(t23);
			t29 = cos(t24);
			t30 = cos(t25);
			t31 = cos(t26);
			t32 = sin(t24);
			t33 = sin(t25);
			t34 = sin(t26);
			t37 = t2.*t3.*t16;
			t38 = t7.*t8.*t11;
			t39 = t7.*t9.*t16;
			t46 = t2.*t3.*t8.*t9;
			t47 = t2.*t3.*t8.*t17;
			t48 = t3.*t8.*t10.*t15;
			t35 = t16.*t18;
			t36 = t16.*t19;
			t40 = t27.^2;
			t41 = t28.^2;
			t42 = t29.^2;
			t43 = t31.^2;
			t51 = -t46;
			t52 = -t48;
			t55 = t20+t39;
			t63 = t29.*t30.*t31.*t32.*t33.*t34.*8.0;
			t44 = -t35;
			t45 = -t36;
			t49 = t42.*2.0;
			t50 = t43.*2.0;
			t58 = t11.*t55;
			t59 = t42.*t43.*4.0;
			t62 = t37+t38+t52;
			t53 = -t49;
			t54 = -t50;
			t56 = t18+t45;
			t57 = t19+t44;
			t60 = t11.*t57;
			t61 = t3.*t10.*t56;
			t65 = t53+t54+t59+t63+1.0;
			t64 = t51+t58+t61;
			t66 = t3.*t10.*t65;
			t67 = t47+t60+t66;
			out1 = reshape([0.0,0.0,0.0,-t62.*(t12.*t14+t4.*t6.*t13)-t67.*(t4.*t14-t6.*t12.*t13)-t5.*t6.*t64,t62.*(t6.*t12-t4.*t13.*t14)+t67.*(t40.*-2.0-t41.*2.0+t40.*t41.*4.0+t27.*t28.*cos(t22).*sin(t21).*sin(t22).*sin(t23).*8.0+1.0)-t5.*t14.*t64,t13.*t64-t4.*t5.*t62+t5.*t12.*t67,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_1(t,in2,in3,in4,in5,in6)
			%InputFVJac_1
			%    OUT1 = InputFVJac_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_2(t,in2,in3,in4,in5,in6)
			%InputFVJac_2
			%    OUT1 = InputFVJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:07
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_3(t,in2,in3,in4,in5,in6)
			%InputFVJac_3
			%    OUT1 = InputFVJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_4(t,in2,in3,in4,in5,in6)
			%InputFVJac_4
			%    OUT1 = InputFVJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:08
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_5(t,in2,in3,in4,in5,in6)
			%InputFVJac_5
			%    OUT1 = InputFVJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:09
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
			function out1 = InputFVJac_6(t,in2,in3,in4,in5,in6)
			%InputFVJac_6
			%    OUT1 = InputFVJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:10
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
		end
		function [ConsJac,ConsCor,ConsGF]=Constraint_Arm(t,q,p,u,s,TransDis_Global,Vel_Global,Cor_Global,Jac_Global,Quat_Global)
		%% Constraint dynamic property calculator (Toolbox Internal Use)
		%% 
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    ConsNum=1;
		    if(ConsNum>0)
		        ConsJac=zeros(ConsNum,numel(q)/2);
		        ConsCor=zeros(ConsNum,1);
		        ConsGF=zeros(ConsNum,1);
		    else
		        ConsJac=zeros(1,numel(q)/2);
		        ConsCor=0;
		        ConsGF=0;
		    end
		    
		    SubConsJac=zeros(1,numel(q)/2);
		    SubConsCor=0;
		    SubConsGF=0;
		    for ConsCount=1:ConsNum
		        
		%SWITCHCASE_
				switch ConsCount
				    case 1
				        SubFrame=[0];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_1(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_1(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_1(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_1(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            CorSubs=Cor_Global(:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
				                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
				            end
				        end
				end
		        
		        ConsJac(ConsCount,:)=SubConsJac;
		        ConsCor(ConsCount,:)=SubConsCor;
		        ConsGF(ConsCount,:)=SubConsGF;
		    end
			function out1 = ConsJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsJac_1
			%    OUT1 = ConsJac_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:10
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			xx3__dt_0_ = in5(1,:);
			zz3__dt_0_ = in5(3,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = sin(q1dev__dt_0_);
			t6 = sin(q1flex__dt_0_);
			t7 = sin(q1sup__dt_0_);
			t9 = q1dev__dt_0_./2.0;
			t10 = q1flex__dt_0_./2.0;
			t11 = q1sup__dt_0_./2.0;
			t12 = xx3__dt_0_./2.0;
			t13 = zz3__dt_0_./2.0;
			t8 = -t4;
			t14 = cos(t9);
			t15 = cos(t10);
			t16 = cos(t11);
			t17 = cos(t12);
			t18 = cos(t13);
			t19 = sin(t9);
			t20 = sin(t10);
			t21 = sin(t11);
			t22 = t7./2.0;
			t23 = sin(t12);
			t24 = sin(t13);
			t25 = (t2.*t4)./2.0;
			t26 = (t3.*t4)./2.0;
			t27 = t14.^2;
			t28 = t15.^2;
			t29 = t16.^2;
			t30 = t16.^3;
			t39 = t16.*t19.*t20;
			t40 = t14.*t15.*t16;
			t41 = t19.*t20.*t21.*2.0;
			t44 = (t14.*t15.*t21)./2.0;
			t45 = (t19.*t20.*t21)./2.0;
			t50 = t4.*t19.*t20.*t21.*(3.0./2.0);
			t31 = t29.^2;
			t32 = t29.*3.0;
			t34 = -t27;
			t35 = -t28;
			t37 = t27.*t29.*2.0;
			t38 = t28.*t29.*2.0;
			t42 = t19.*t20.*t30;
			t43 = -t40;
			t46 = t14.*t15.*t30.*2.0;
			t47 = t14.*t15.*t30.*3.0;
			t49 = t4.*t44;
			t51 = -t50;
			t52 = t29.*t41;
			t53 = t19.*t20.*t21.*t29.*-2.0;
			t33 = t31.*2.0;
			t36 = -t32;
			t48 = -t42;
			t54 = t22+t39+t44+t48+t49;
			t55 = t8+t25+t26+t33+t43+t45+t47+t51;
			t56 = t31+t34+t35+t36+t37+t38+t41+t46+t53+2.0;
			t57 = 1.0./t56;
			out1 = [t14.*t16.*t20.*(-1.0./2.0)+(t15.*t19.*t21)./2.0-(t17.*t24.*t54.*t57)./4.0+(t18.*t23.*t55.*t57)./4.0,t15.*t16.*t19.*(-1.0./2.0)+(t14.*t20.*t21)./2.0-(t18.*t23.*t54.*t57)./4.0+(t17.*t24.*t55.*t57)./4.0,t40.*(-1.0./2.0)+t45-(t17.*t24.*t57.*(t5./2.0+t6.*t22+t15.*t16.*t19.*2.0+t14.*t20.*t21.*(3.0./2.0)-t15.*t19.*t30+(t4.*t14.*t20.*t21)./2.0))./4.0-(t18.*t23.*t57.*(t6./2.0+t5.*t22+t14.*t16.*t20.*2.0+t15.*t19.*t21.*(3.0./2.0)-t14.*t20.*t30+(t4.*t15.*t19.*t21)./2.0))./4.0];
			end
			function out1 = ConsCor_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsCor_1
			%    OUT1 = ConsCor_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:12
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(4,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(5,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(6,:);
			xx3__dt_0_ = in5(1,:);
			zz3__dt_0_ = in5(3,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = sin(q1dev__dt_0_);
			t6 = sin(q1flex__dt_0_);
			t7 = sin(q1sup__dt_0_);
			t12 = q1dev__dt_0_./2.0;
			t13 = q1flex__dt_0_./2.0;
			t14 = q1sup__dt_0_./2.0;
			t15 = xx3__dt_0_./2.0;
			t16 = zz3__dt_0_./2.0;
			t8 = q1dev__dt_1_.*t4;
			t9 = q1flex__dt_1_.*t4;
			t10 = -t4;
			t11 = -t7;
			t17 = cos(t12);
			t18 = cos(t13);
			t19 = cos(t14);
			t20 = t4./2.0;
			t21 = cos(t15);
			t22 = cos(t16);
			t23 = sin(t12);
			t24 = t5./2.0;
			t25 = sin(t13);
			t26 = t6./2.0;
			t27 = sin(t14);
			t28 = t7./2.0;
			t29 = sin(t15);
			t30 = sin(t16);
			t31 = q1dev__dt_1_.*t28;
			t32 = q1sup__dt_1_.*t24;
			t33 = q1flex__dt_1_.*t28;
			t34 = q1sup__dt_1_.*t26;
			t35 = t2.*t20;
			t36 = t3.*t20;
			t37 = t2.*t28;
			t38 = t5.*t20;
			t39 = t3.*t28;
			t40 = t6.*t20;
			t41 = t7.*t24;
			t42 = t7.*t26;
			t43 = t17.^2;
			t44 = t18.^2;
			t45 = t19.^2;
			t46 = t19.^3;
			t48 = t27.^2;
			t49 = (t2.*t8)./2.0;
			t50 = (t3.*t8)./2.0;
			t51 = (t2.*t9)./2.0;
			t52 = (t3.*t9)./2.0;
			t68 = t17.*t23;
			t69 = t18.*t25;
			t70 = t19.*t27.*3.0;
			t78 = t17.*t25.*t27;
			t79 = t18.*t23.*t27;
			t80 = t19.*t23.*t25;
			t81 = t17.*t18.*t19;
			t82 = t17.*t19.*t25.*2.0;
			t83 = t18.*t19.*t23.*2.0;
			t84 = t23.*t25.*t27.*2.0;
			t100 = (t17.*t18.*t27)./2.0;
			t101 = (t17.*t19.*t25)./2.0;
			t102 = (t18.*t19.*t23)./2.0;
			t103 = (t17.*t18.*t27)./4.0;
			t104 = (t17.*t19.*t25)./4.0;
			t105 = (t18.*t19.*t23)./4.0;
			t106 = t17.*t18.*t27.*(3.0./4.0);
			t112 = (t23.*t25.*t27)./2.0;
			t113 = (t23.*t25.*t27)./4.0;
			t114 = t23.*t25.*t27.*(3.0./4.0);
			t156 = t7.*t23.*t25.*t27.*(3.0./2.0);
			t165 = q1dev__dt_1_.*t23.*t25.*t27.*(-1.0./2.0);
			t166 = q1flex__dt_1_.*t23.*t25.*t27.*(-1.0./2.0);
			t172 = t17.*t18.*t20.*t27;
			t176 = t17.*t18.*t27.*t28;
			t185 = t4.*t23.*t25.*t27.*(3.0./2.0);
			t192 = t8.*t23.*t25.*t27.*(3.0./2.0);
			t193 = t9.*t23.*t25.*t27.*(3.0./2.0);
			t197 = t7.*t17.*t18.*t27.*(-1.0./2.0);
			t47 = t45.^2;
			t53 = t7.*t32;
			t54 = t7.*t34;
			t55 = -t49;
			t56 = -t50;
			t57 = -t51;
			t58 = -t52;
			t59 = t45.*3.0;
			t63 = -t43;
			t64 = -t44;
			t71 = -t68;
			t72 = -t69;
			t73 = -t70;
			t74 = t27.*t46.*2.0;
			t75 = t27.*t46.*4.0;
			t76 = t43.*t45.*2.0;
			t77 = t44.*t45.*2.0;
			t85 = q1dev__dt_1_.*t81;
			t86 = q1flex__dt_1_.*t81;
			t87 = q1sup__dt_1_.*t81;
			t88 = q1dev__dt_1_.*t80;
			t89 = q1flex__dt_1_.*t80;
			t90 = q1sup__dt_1_.*t80;
			t91 = -t78;
			t92 = -t79;
			t93 = -t80;
			t94 = t17.*t25.*t46;
			t95 = t18.*t23.*t46;
			t96 = t23.*t25.*t46;
			t97 = q1sup__dt_1_.*t82;
			t98 = q1sup__dt_1_.*t83;
			t99 = -t81;
			t107 = t78./4.0;
			t108 = t78.*(3.0./2.0);
			t109 = t79./4.0;
			t110 = t79.*(3.0./2.0);
			t111 = t80./4.0;
			t115 = t17.*t18.*t46.*2.0;
			t116 = t17.*t18.*t46.*3.0;
			t117 = t45.*t68.*2.0;
			t118 = t45.*t69.*2.0;
			t119 = t19.*t27.*t43.*2.0;
			t120 = t19.*t27.*t44.*2.0;
			t127 = t81./4.0;
			t128 = -t100;
			t129 = -t101;
			t130 = -t102;
			t134 = -t112;
			t138 = q1dev__dt_1_.*t100;
			t139 = q1dev__dt_1_.*t101;
			t140 = q1dev__dt_1_.*t102;
			t141 = q1flex__dt_1_.*t100;
			t142 = q1flex__dt_1_.*t101;
			t143 = q1flex__dt_1_.*t102;
			t144 = q1sup__dt_1_.*t106;
			t151 = q1dev__dt_1_.*t112;
			t152 = q1flex__dt_1_.*t112;
			t153 = q1sup__dt_1_.*t114;
			t157 = (t17.*t18.*t46)./2.0;
			t163 = q1dev__dt_1_.*t78.*(-1.0./4.0);
			t164 = q1flex__dt_1_.*t79.*(-1.0./4.0);
			t167 = q1dev__dt_1_.*t17.*t18.*t46.*-3.0;
			t168 = q1flex__dt_1_.*t17.*t18.*t46.*-3.0;
			t173 = t4.*t103;
			t177 = t20.*t78;
			t178 = t20.*t79;
			t181 = t4.*t78.*(3.0./4.0);
			t182 = t4.*t79.*(3.0./4.0);
			t183 = t4.*t80.*(3.0./4.0);
			t184 = t4.*t113;
			t186 = -t156;
			t187 = t8.*t100;
			t188 = t9.*t100;
			t198 = -t185;
			t199 = t45.*t78;
			t200 = t45.*t79;
			t202 = t48.*t80.*2.0;
			t203 = t45.*t84;
			t205 = t23.*t25.*t27.*t45.*-2.0;
			t206 = t23.*t25.*t27.*t45.*(3.0./2.0);
			t207 = t17.*t18.*t27.*t45.*(9.0./2.0);
			t60 = t47.*2.0;
			t65 = -t59;
			t66 = q1dev__dt_1_.*t47.*-2.0;
			t67 = q1flex__dt_1_.*t47.*-2.0;
			t121 = -t87;
			t122 = -t90;
			t123 = q1sup__dt_1_.*t94;
			t124 = q1sup__dt_1_.*t95;
			t125 = q1dev__dt_1_.*t96;
			t126 = q1flex__dt_1_.*t96;
			t131 = -t107;
			t132 = -t109;
			t133 = -t111;
			t135 = -t94;
			t136 = -t95;
			t137 = -t96;
			t145 = q1dev__dt_1_.*t107;
			t146 = q1dev__dt_1_.*t109;
			t147 = q1flex__dt_1_.*t107;
			t148 = q1flex__dt_1_.*t109;
			t149 = q1sup__dt_1_.*t108;
			t150 = q1sup__dt_1_.*t110;
			t154 = q1dev__dt_1_.*t116;
			t155 = q1flex__dt_1_.*t116;
			t158 = t94./2.0;
			t159 = t95./2.0;
			t160 = t94.*(3.0./2.0);
			t161 = t95.*(3.0./2.0);
			t162 = t96./2.0;
			t169 = t4.*t127;
			t179 = t4.*t107;
			t180 = t4.*t109;
			t189 = q1sup__dt_1_.*t173;
			t190 = q1sup__dt_1_.*t177;
			t191 = q1sup__dt_1_.*t178;
			t194 = q1sup__dt_1_.*t184;
			t195 = q1sup__dt_1_.*t157;
			t201 = t17.*t18.*t27.*t59;
			t204 = -t202;
			t208 = t103+t111;
			t209 = t104+t109;
			t210 = t105+t107;
			t211 = t113+t127;
			t219 = t71+t91+t95+t117+t199;
			t220 = t72+t92+t94+t118+t200;
			t61 = q1dev__dt_1_.*t60;
			t62 = q1flex__dt_1_.*t60;
			t170 = -t123;
			t171 = -t124;
			t174 = -t125;
			t175 = -t126;
			t196 = q1sup__dt_1_.*t162;
			t212 = t109+t129+t158+t180;
			t213 = t107+t130+t159+t179;
			t214 = t28+t80+t100+t137+t172;
			t215 = t38+t130+t131+t161+t181;
			t216 = t40+t129+t132+t160+t182;
			t217 = t24+t42+t83+t108+t136+t177;
			t218 = t26+t41+t82+t110+t135+t178;
			t221 = t10+t35+t36+t60+t99+t112+t116+t198;
			t222 = t20+t127+t134+t169+t197+t206;
			t223 = t47+t63+t64+t65+t76+t77+t84+t115+t205+2.0;
			t226 = t11+t37+t39+t75+t128+t133+t183+t186+t207;
			t227 = t73+t74+t93+t96+t119+t120+t201+t204;
			t224 = 1.0./t223;
			t232 = t8+t33+t34+t53+t55+t56+t66+t85+t89+t97+t141+t150+t165+t167+t170+t175+t188+t191+t192;
			t233 = t9+t31+t32+t54+t57+t58+t67+t86+t88+t98+t138+t149+t166+t168+t171+t174+t187+t190+t193;
			t225 = t224.^2;
			t228 = (t21.*t22.*t214.*t224)./8.0;
			t229 = (t29.*t30.*t214.*t224)./8.0;
			t230 = (t21.*t22.*t221.*t224)./8.0;
			t231 = (t29.*t30.*t221.*t224)./8.0;
			t234 = (t21.*t22.*t224.*t232)./8.0;
			t235 = (t21.*t22.*t224.*t233)./8.0;
			t236 = (t29.*t30.*t224.*t232)./8.0;
			t237 = (t29.*t30.*t224.*t233)./8.0;
			t238 = -t236;
			t239 = -t237;
			t240 = t229+t230;
			t241 = t228+t231;
			t242 = t234+t239;
			t243 = t235+t238;
			et1 = -q1dev__dt_1_.*t208+q1flex__dt_1_.*t211-q1sup__dt_1_.*t210-q1dev__dt_1_.*(t208+(t21.*t30.*t212.*t224)./4.0-(t22.*t29.*t215.*t224)./4.0-(t21.*t30.*t214.*t219.*t225)./4.0+(t22.*t29.*t219.*t221.*t225)./4.0)+q1flex__dt_1_.*(t211-(t21.*t30.*t213.*t224)./4.0+(t22.*t29.*t216.*t224)./4.0+(t21.*t30.*t214.*t220.*t225)./4.0-(t22.*t29.*t220.*t221.*t225)./4.0)+q1sup__dt_1_.*(-t210+(t21.*t30.*t222.*t224)./4.0+(t22.*t29.*t224.*t226)./4.0+(t21.*t30.*t214.*t225.*t227)./4.0-(t22.*t29.*t221.*t225.*t227)./4.0)-t214.*t224.*t242+t221.*t224.*t243-t224.*t232.*t241+t224.*t233.*t240;
			et2 = t21.*t30.*t224.*(t121+t143+t146+t147+t153+t194+t195+q1dev__dt_1_.*t158-q1flex__dt_1_.*t95.*(3.0./2.0)-(q1sup__dt_1_.*t2)./2.0-(t5.*t9)./2.0-t9.*t78.*(3.0./4.0)+t8.*t109-(q1dev__dt_1_.*t17.*t19.*t25)./2.0).*(-1.0./4.0)+(t22.*t29.*t224.*(t122+t142+t144+t163+t164+t189+t196+q1dev__dt_1_.*t161-(q1flex__dt_1_.*t94)./2.0+q1sup__dt_1_.*t37+t8.*t24+t8.*t78.*(3.0./4.0)-(t9.*t79)./4.0-(q1dev__dt_1_.*t18.*t19.*t23)./2.0))./4.0+(t22.*t29.*t219.*t225.*t232)./4.0+(t21.*t30.*t219.*t225.*t233)./4.0;
			et3 = q1dev__dt_1_.*t211-q1flex__dt_1_.*t208-q1sup__dt_1_.*t209+q1dev__dt_1_.*(t211-(t22.*t29.*t212.*t224)./4.0+(t21.*t30.*t215.*t224)./4.0+(t22.*t29.*t214.*t219.*t225)./4.0-(t21.*t30.*t219.*t221.*t225)./4.0)-q1flex__dt_1_.*(t208+(t22.*t29.*t213.*t224)./4.0-(t21.*t30.*t216.*t224)./4.0-(t22.*t29.*t214.*t220.*t225)./4.0+(t21.*t30.*t220.*t221.*t225)./4.0)+q1sup__dt_1_.*(-t209+(t22.*t29.*t222.*t224)./4.0+(t21.*t30.*t224.*t226)./4.0+(t22.*t29.*t214.*t225.*t227)./4.0-(t21.*t30.*t221.*t225.*t227)./4.0)-t214.*t224.*t243+t221.*t224.*t242+t224.*t232.*t240-t224.*t233.*t241;
			et4 = t22.*t29.*t224.*(t121+t139+t146+t147+t153+t194+t195-q1dev__dt_1_.*t94.*(3.0./2.0)+q1flex__dt_1_.*t159-(q1sup__dt_1_.*t3)./2.0-(t6.*t8)./2.0-t8.*t79.*(3.0./4.0)+t9.*t107-(q1flex__dt_1_.*t18.*t19.*t23)./2.0).*(-1.0./4.0)+(t21.*t30.*t224.*(t122+t140+t144+t163+t164+t189+t196-(q1dev__dt_1_.*t95)./2.0+q1flex__dt_1_.*t160+q1sup__dt_1_.*t39+t9.*t26-(t8.*t78)./4.0+t9.*t79.*(3.0./4.0)-(q1flex__dt_1_.*t17.*t19.*t25)./2.0))./4.0+(t22.*t29.*t220.*t225.*t232)./4.0+(t21.*t30.*t220.*t225.*t233)./4.0;
			et5 = q1dev__dt_1_.*t210+q1flex__dt_1_.*t209+q1sup__dt_1_.*t208-q1dev__dt_1_.*(-t210+(t22.*t29.*t224.*(t37+t93+t106+t162+t173))./4.0-(t21.*t30.*t224.*(t2.*(-1.0./2.0)+t99+t114+t157+t184))./4.0+(t21.*t30.*t217.*t219.*t225)./4.0+(t22.*t29.*t218.*t219.*t225)./4.0)-q1flex__dt_1_.*(-t209+(t21.*t30.*t224.*(t39+t93+t106+t162+t173))./4.0-(t22.*t29.*t224.*(t3.*(-1.0./2.0)+t99+t114+t157+t184))./4.0+(t21.*t30.*t217.*t220.*t225)./4.0+(t22.*t29.*t218.*t220.*t225)./4.0);
			et6 = -q1sup__dt_1_.*(-t103+t133+(t22.*t29.*t224.*(t38+t91-(t7.*t79)./2.0+t4.*t105+t45.*t108+t18.*t19.*t23.*(3.0./4.0)))./4.0+(t21.*t30.*t224.*(t40+t92-(t7.*t78)./2.0+t4.*t104+t45.*t110+t17.*t19.*t25.*(3.0./4.0)))./4.0+(t21.*t30.*t217.*t225.*t227)./4.0+(t22.*t29.*t218.*t225.*t227)./4.0)+t224.*t232.*((t21.*t22.*t217.*t224)./8.0-(t29.*t30.*t218.*t224)./8.0)+t224.*t233.*((t21.*t22.*t218.*t224)./8.0-(t29.*t30.*t217.*t224)./8.0)+t217.*t224.*t242+t218.*t224.*t243;
			et7 = t22.*t29.*t224.*(t9./2.0+t86./4.0-t88./4.0+t166+q1dev__dt_1_.*t11+q1dev__dt_1_.*t75+q1dev__dt_1_.*t207+q1flex__dt_1_.*t197+q1flex__dt_1_.*t206+q1sup__dt_1_.*t38+q1sup__dt_1_.*t91+t2.*t31+t3.*t31+t8.*t80.*(3.0./4.0)+t9.*t127+t45.*t149-(q1sup__dt_1_.*t7.*t79)./2.0+q1sup__dt_1_.*t4.*t105-(q1dev__dt_1_.*t17.*t18.*t27)./2.0+q1sup__dt_1_.*t18.*t19.*t23.*(3.0./4.0)-q1dev__dt_1_.*t7.*t23.*t25.*t27.*(3.0./2.0)).*(-1.0./4.0);
			et8 = t21.*t30.*t224.*(t8./2.0+t85./4.0-t89./4.0+t165+q1dev__dt_1_.*t197+q1dev__dt_1_.*t206+q1flex__dt_1_.*t11+q1flex__dt_1_.*t75+q1flex__dt_1_.*t207+q1sup__dt_1_.*t40+q1sup__dt_1_.*t92+t2.*t33+t3.*t33+t9.*t80.*(3.0./4.0)+t8.*t127+t45.*t150-(q1sup__dt_1_.*t7.*t78)./2.0+q1sup__dt_1_.*t4.*t104-(q1flex__dt_1_.*t17.*t18.*t27)./2.0+q1sup__dt_1_.*t17.*t19.*t25.*(3.0./4.0)-q1flex__dt_1_.*t7.*t23.*t25.*t27.*(3.0./2.0)).*(-1.0./4.0)-(t22.*t29.*t225.*t227.*t232)./4.0-(t21.*t30.*t225.*t227.*t233)./4.0;
			out1 = q1sup__dt_1_.*(et5+et6+et7+et8)-q1dev__dt_1_.*(et1+et2)-q1flex__dt_1_.*(et3+et4);
			end
			function out1 = ConsGF_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsGF_1
			%    OUT1 = ConsGF_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:14
			c = in3(3,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(4,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(5,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(6,:);
			xx3__dt_0_ = in5(1,:);
			zz3__dt_0_ = in5(3,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = sin(q1dev__dt_0_);
			t6 = sin(q1flex__dt_0_);
			t7 = sin(q1sup__dt_0_);
			t8 = q1dev__dt_0_./2.0;
			t9 = q1flex__dt_0_./2.0;
			t10 = q1sup__dt_0_./2.0;
			t11 = xx3__dt_0_./2.0;
			t12 = zz3__dt_0_./2.0;
			t13 = cos(t8);
			t14 = cos(t9);
			t15 = cos(t10);
			t16 = sin(t8);
			t17 = sin(t9);
			t18 = sin(t10);
			t19 = sin(t11);
			t20 = sin(t12);
			t21 = t13.^2;
			t22 = t14.^2;
			t23 = t15.^2;
			t24 = t15.^3;
			t32 = t16.*t17.*t18.*2.0;
			t25 = t23.^2;
			t26 = t23.*3.0;
			t27 = -t21;
			t28 = -t22;
			t30 = t21.*t23.*2.0;
			t31 = t22.*t23.*2.0;
			t33 = t13.*t14.*t24.*2.0;
			t34 = t23.*t32;
			t35 = t16.*t17.*t18.*t23.*-2.0;
			t29 = -t26;
			t36 = t25+t27+t28+t29+t30+t31+t32+t33+t35+2.0;
			t37 = 1.0./t36;
			et1 = q1dev__dt_1_.*((t13.*t15.*t17)./2.0-(t14.*t16.*t18)./2.0)+q1flex__dt_1_.*((t14.*t15.*t16)./2.0-(t13.*t17.*t18)./2.0)+q1sup__dt_1_.*((t13.*t14.*t15)./2.0-(t16.*t17.*t18)./2.0)+(t19.*t37.*cos(t12).*(q1dev__dt_1_.*t4-q1dev__dt_1_.*t25.*2.0+(q1flex__dt_1_.*t7)./2.0+(q1sup__dt_1_.*t6)./2.0-(q1dev__dt_1_.*t2.*t4)./2.0-(q1dev__dt_1_.*t3.*t4)./2.0+(q1sup__dt_1_.*t5.*t7)./2.0+q1dev__dt_1_.*t13.*t14.*t15-q1dev__dt_1_.*t13.*t14.*t24.*3.0-(q1dev__dt_1_.*t16.*t17.*t18)./2.0+(q1flex__dt_1_.*t13.*t14.*t18)./2.0+q1flex__dt_1_.*t15.*t16.*t17-q1flex__dt_1_.*t16.*t17.*t24+q1sup__dt_1_.*t13.*t15.*t17.*2.0+q1sup__dt_1_.*t14.*t16.*t18.*(3.0./2.0)-q1sup__dt_1_.*t13.*t17.*t24+q1dev__dt_1_.*t4.*t16.*t17.*t18.*(3.0./2.0)+(q1flex__dt_1_.*t4.*t13.*t14.*t18)./2.0+(q1sup__dt_1_.*t4.*t14.*t16.*t18)./2.0))./4.0;
			et2 = (t20.*t37.*cos(t11).*((q1dev__dt_1_.*t7)./2.0+q1flex__dt_1_.*t4-q1flex__dt_1_.*t25.*2.0+(q1sup__dt_1_.*t5)./2.0-(q1flex__dt_1_.*t2.*t4)./2.0-(q1flex__dt_1_.*t3.*t4)./2.0+(q1sup__dt_1_.*t6.*t7)./2.0+(q1dev__dt_1_.*t13.*t14.*t18)./2.0+q1dev__dt_1_.*t15.*t16.*t17-q1dev__dt_1_.*t16.*t17.*t24+q1flex__dt_1_.*t13.*t14.*t15-q1flex__dt_1_.*t13.*t14.*t24.*3.0-(q1flex__dt_1_.*t16.*t17.*t18)./2.0+q1sup__dt_1_.*t14.*t15.*t16.*2.0+q1sup__dt_1_.*t13.*t17.*t18.*(3.0./2.0)-q1sup__dt_1_.*t14.*t16.*t24+(q1dev__dt_1_.*t4.*t13.*t14.*t18)./2.0+q1flex__dt_1_.*t4.*t16.*t17.*t18.*(3.0./2.0)+(q1sup__dt_1_.*t4.*t13.*t17.*t18)./2.0))./4.0;
			out1 = -c.*(t19.*t20.*(-1.0./2.0)+t13.*t14.*t18+t15.*t16.*t17)-c.*(et1+et2);
			end
			function out1 = ConsFVJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsFVJac_1
			%    OUT1 = ConsFVJac_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 20:12:14
			out1 = 0.0;
			end
		end
	end

end

