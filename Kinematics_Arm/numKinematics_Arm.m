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

