function [flow,flowTraj,q,p,s]=Jump_TAWE(flow,t,q,p,u,s,l)
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

	function [TF,Vel,Cor,Jac,TransDis,Quat,MM,GFI,GFF,GFU,JacU,JacCons,CorCons,GFCons,CenMat,CorMatLeft,CorMatRight]=System_TAWE(t,q,p,u,s)
	%% Author: Jiamin Wang; Updated: 2021-12-15;
	    
	    [TF,Vel,Cor,Jac,TransDis,Quat]=numKinematics_TAWE(t,q,p,u,s);
	    [MM,GFI,CenMat,CorMatLeft,CorMatRight]=Inertial_TAWE(t,q,p,u,s,TF,Vel,Cor,Jac,true);
	    [GFF]=Force_TAWE(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
	    [GFU,JacU]=Input_TAWE(t,q,p,u,s,TF,TransDis,Vel,Jac,Quat);
	    [JacCons,CorCons,GFCons]=Constraint_TAWE(t,q,p,u,s,TransDis,Vel,Cor,Jac,Quat);
		function [frameTF,frameVel,frameCorAcc,frameJacobian,frameTransDis,frameRotQuat]=numKinematics_TAWE(t,q,p,u,s)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    LinkTF=linkTF(t,q,p,u,s);
		    LinkVel=linkVel(t,q,p,u,s);
		    LinkCorAcc=linkCorAcc(t,q,p,u,s);
		    LinkJacobian=linkJacobian(t,q,p,u,s);
		    LinkRotQuat=linkRotQuat(t,q,p,u,s);
		    FrameNum=26;
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
				    case 12
				        framePath=[1   2  12];
				    case 13
				        framePath=[1   2  12  13];
				    case 14
				        framePath=[1   2  12  13  14];
				    case 15
				        framePath=[1   2  12  13  14  15];
				    case 16
				        framePath=[1   2  12  13  14  15  16];
				    case 17
				        framePath=[1   2  12  13  14  15  16  17];
				    case 18
				        framePath=[1   2  12  13  14  15  16  17  18];
				    case 19
				        framePath=[1   2  12  13  14  15  16  17  18  19];
				    case 20
				        framePath=[1   2  12  13  14  15  16  17  18  19  20];
				    case 21
				        framePath=[1   2  12  13  14  21];
				    case 22
				        framePath=[1   2  12  13  14  15  22];
				    case 23
				        framePath=[1   2  12  13  14  15  16  23];
				    case 24
				        framePath=[1   2  12  13  14  15  16  17  24];
				    case 25
				        framePath=[1   2  12  13  14  15  16  17  18  25];
				    case 26
				        framePath=[1   2  12  13  14  15  16  17  18  19  26];
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
			%    14-Apr-2023 16:19:06
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l11x = in3(7,:);
			l11y = in3(8,:);
			l12x = in3(10,:);
			l11z = in3(9,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			l21x = in3(16,:);
			l21y = in3(17,:);
			l21z = in3(18,:);
			lm13x = in3(42,:);
			lm13y = in3(43,:);
			lm13z = in3(44,:);
			lm22x = in3(45,:);
			lm22y = in3(46,:);
			lm23x = in3(48,:);
			lm22z = in3(47,:);
			lm23y = in3(49,:);
			lm24x = in3(51,:);
			lm23z = in3(50,:);
			lm24y = in3(52,:);
			lm25x = in3(54,:);
			lm24z = in3(53,:);
			lm25y = in3(55,:);
			lm26x = in3(57,:);
			lm25z = in3(56,:);
			lm26y = in3(58,:);
			lm27x = in3(60,:);
			lm26z = in3(59,:);
			lm27y = in3(61,:);
			lm27z = in3(62,:);
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			q2w1__dt_0_ = in2(4,:);
			q2w2__dt_0_ = in2(5,:);
			q2w3__dt_0_ = in2(6,:);
			q2w4__dt_0_ = in2(7,:);
			q2w5__dt_0_ = in2(8,:);
			q2w6__dt_0_ = in2(9,:);
			rTetra = in3(6,:);
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU2x = in3(30,:);
			rf1IMU1z = in3(26,:);
			rf1IMU2y = in3(31,:);
			rf1IMU2z = in3(32,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf1IMU1x);
			t6 = cos(rf1IMU1y);
			t7 = cos(rf1IMU2x);
			t8 = cos(rf1IMU1z);
			t9 = cos(rf1IMU2y);
			t10 = cos(rf1IMU2z);
			t11 = cos(rf2Medx);
			t12 = cos(rf2Medy);
			t13 = cos(rf2Medz);
			t14 = cos(rf1Devx);
			t15 = cos(rf1Devy);
			t16 = cos(rf1Devz);
			t17 = sin(q1dev__dt_0_);
			t18 = sin(q1flex__dt_0_);
			t19 = sin(q1sup__dt_0_);
			t20 = sin(rf1IMU1x);
			t21 = sin(rf1IMU1y);
			t22 = sin(rf1IMU2x);
			t23 = sin(rf1IMU1z);
			t24 = sin(rf1IMU2y);
			t25 = sin(rf1IMU2z);
			t26 = sin(rf2Medx);
			t27 = sin(rf2Medy);
			t28 = sin(rf2Medz);
			t29 = sin(rf1Devx);
			t30 = sin(rf1Devy);
			t31 = sin(rf1Devz);
			t32 = q2w1__dt_0_+riw1;
			t33 = q2w2__dt_0_+riw2;
			t34 = q2w3__dt_0_+riw3;
			t35 = q2w4__dt_0_+riw4;
			t36 = q2w5__dt_0_+riw5;
			t37 = q2w6__dt_0_+riw6;
			t58 = sqrt(3.0);
			t59 = sqrt(6.0);
			t60 = q1dev__dt_0_./2.0;
			t61 = q1flex__dt_0_./2.0;
			t62 = q1sup__dt_0_./2.0;
			t63 = rf1IMU1x./2.0;
			t64 = rf1IMU1y./2.0;
			t65 = rf1IMU2x./2.0;
			t66 = rf1IMU1z./2.0;
			t67 = rf1IMU2y./2.0;
			t68 = rf1IMU2z./2.0;
			t69 = rf2Medx./2.0;
			t70 = rf2Medy./2.0;
			t71 = rf2Medz./2.0;
			t72 = rf1Devx./2.0;
			t73 = rf1Devy./2.0;
			t74 = rf1Devz./2.0;
			t38 = cos(t32);
			t39 = cos(t33);
			t40 = cos(t34);
			t41 = cos(t35);
			t42 = cos(t36);
			t43 = cos(t37);
			t44 = sin(t32);
			t45 = sin(t33);
			t46 = sin(t34);
			t47 = sin(t35);
			t48 = sin(t36);
			t49 = sin(t37);
			t50 = t7.*t9;
			t51 = t9.*t10;
			t52 = t9.*t22;
			t53 = t7.*t25;
			t54 = t10.*t22;
			t55 = t9.*t25;
			t56 = t22.*t25;
			t57 = -t24;
			t75 = cos(t60);
			t76 = cos(t61);
			t77 = cos(t63);
			t78 = cos(t65);
			t79 = cos(t66);
			t80 = cos(t67);
			t81 = cos(t68);
			t82 = cos(t69);
			t83 = cos(t71);
			t84 = cos(t72);
			t85 = cos(t74);
			t86 = sin(t65);
			t87 = sin(t67);
			t88 = sin(t68);
			t93 = t7.*t10.*t24;
			t104 = (rTetra.*t58)./3.0;
			t105 = (rTetra.*t59)./6.0;
			t89 = t24.*t53;
			t90 = t24.*t54;
			t91 = -t53;
			t92 = -t54;
			t94 = t75.^2;
			t95 = t76.^2;
			t96 = t77.^2;
			t97 = t78.^2;
			t98 = t79.^2;
			t99 = t81.^2;
			t100 = t82.^2;
			t101 = t83.^2;
			t102 = t84.^2;
			t103 = t85.^2;
			t108 = -t104;
			t109 = -t105;
			t112 = t56+t93;
			t116 = t78.*t80.*t81.*t86.*t87.*t88.*8.0;
			t106 = t97.*2.0;
			t107 = t99.*2.0;
			t113 = t90+t91;
			t114 = t89+t92;
			t115 = t97.*t99.*4.0;
			t110 = -t106;
			t111 = -t107;
			t117 = t110+t111+t115+t116+1.0;
			mt1 = [1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,t6.*t8,t6.*t23,-t21,0.0,-t5.*t23+t8.*t20.*t21,t96.*-2.0-t98.*2.0+t96.*t98.*4.0+t77.*t79.*cos(t64).*sin(t63).*sin(t64).*sin(t66).*8.0+1.0,t6.*t20,0.0,t20.*t23+t5.*t8.*t21,-t8.*t20+t5.*t21.*t23,t5.*t6,0.0,0.0,0.0,0.0,1.0,t15.*t16,t15.*t31,-t30,0.0,-t14.*t31+t16.*t29.*t30,t102.*-2.0-t103.*2.0+t102.*t103.*4.0+t84.*t85.*cos(t73).*sin(t72).*sin(t73).*sin(t74).*8.0+1.0,t15.*t29,0.0,t29.*t31+t14.*t16.*t30,-t16.*t29+t14.*t30.*t31,t14.*t15,0.0,l11x,l11y,l11z,1.0,t2.*t4,t4.*t17,-t19,0.0,-t3.*t17+t2.*t18.*t19];
			mt2 = [t94.*-2.0-t95.*2.0+t94.*t95.*4.0+t75.*t76.*cos(t62).*sin(t60).*sin(t61).*sin(t62).*8.0+1.0,t4.*t18,0.0,t17.*t18+t2.*t3.*t19,-t2.*t18+t3.*t17.*t19,t3.*t4,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,l12x,l12y,l12z,1.0,t51,t55,t57,0.0,t113,t117,t52,0.0,t112,t114,t50,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm13x,lm13y,lm13z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,-rTetra,t109,t108,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,rTetra,t109,t108,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,t109,rTetra.*t58.*(2.0./3.0),1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,(rTetra.*t59)./2.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0];
			mt3 = [t12.*t13,t12.*t28,-t27,0.0,-t11.*t28+t13.*t26.*t27,t100.*-2.0-t101.*2.0+t100.*t101.*4.0+t82.*t83.*cos(t70).*sin(t69).*sin(t70).*sin(t71).*8.0+1.0,t12.*t26,0.0,t26.*t28+t11.*t13.*t27,-t13.*t26+t11.*t27.*t28,t11.*t12,0.0,l21x,l21y,l21z,1.0,t38,t44,0.0,0.0,-t44,t38,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,t39,t45,0.0,0.0,-t45,t39,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,t40,t46,0.0,0.0,-t46,t40,0.0,0.0,l21,0.0,1.0,1.0,0.0,0.0,0.0,0.0,t41,t47,0.0,0.0,-t47,t41,0.0,0.0,l22,0.0,1.0,t42,0.0,-t48,0.0,0.0,1.0,0.0,0.0,t48,0.0,t42,0.0,0.0,l23,0.0,1.0,t43,t49,0.0,0.0,-t49,t43,0.0,0.0,0.0,0.0,1.0,0.0,l24,0.0,l25,1.0,t51,t113,t112,0.0,t55,t117,t114,0.0];
			mt4 = [t57,t52,t50,0.0,l13x,l13y,l13z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm22x,lm22y,lm22z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm23x,lm23y,lm23z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm24x,lm24y,lm24z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm25x,lm25y,lm25z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm26x,lm26y,lm26z,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,lm27x,lm27y,lm27z,1.0];
			out1 = reshape([mt1,mt2,mt3,mt4],4,104);
			end
			function out1 = linkVel(t,in2,in3,in4,in5)
			%linkVel
			%    OUT1 = linkVel(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:19:05
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_1_ = in2(18,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = sin(q1dev__dt_0_);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-q1sup__dt_1_.*t4+q1flex__dt_1_.*t2.*t3,q1sup__dt_1_.*t2+q1flex__dt_1_.*t3.*t4,q1dev__dt_1_-q1flex__dt_1_.*sin(q1sup__dt_0_),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,q2w1__dt_1_,0.0,0.0,0.0,q2w2__dt_1_,0.0,0.0,0.0,0.0,0.0,q2w3__dt_1_,0.0,0.0,0.0,0.0,0.0,q2w4__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,q2w5__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,q2w6__dt_1_,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,26]);
			end
			function out1 = linkRotQuat(t,in2,in3,in4,in5)
			%linkRotQuat
			%    OUT1 = linkRotQuat(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:19:05
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			q2w1__dt_0_ = in2(4,:);
			q2w2__dt_0_ = in2(5,:);
			q2w3__dt_0_ = in2(6,:);
			q2w4__dt_0_ = in2(7,:);
			q2w5__dt_0_ = in2(8,:);
			q2w6__dt_0_ = in2(9,:);
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU2x = in3(30,:);
			rf1IMU1z = in3(26,:);
			rf1IMU2y = in3(31,:);
			rf1IMU2z = in3(32,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = q1dev__dt_0_./2.0;
			t3 = q1flex__dt_0_./2.0;
			t4 = q1sup__dt_0_./2.0;
			t5 = q2w1__dt_0_./2.0;
			t6 = q2w2__dt_0_./2.0;
			t7 = q2w3__dt_0_./2.0;
			t8 = q2w4__dt_0_./2.0;
			t9 = q2w5__dt_0_./2.0;
			t10 = q2w6__dt_0_./2.0;
			t11 = rf1IMU1x./2.0;
			t12 = rf1IMU1y./2.0;
			t13 = rf1IMU2x./2.0;
			t14 = rf1IMU1z./2.0;
			t15 = rf1IMU2y./2.0;
			t16 = rf1IMU2z./2.0;
			t17 = rf2Medx./2.0;
			t18 = rf2Medy./2.0;
			t19 = rf2Medz./2.0;
			t20 = rf1Devx./2.0;
			t21 = rf1Devy./2.0;
			t22 = rf1Devz./2.0;
			t23 = riw1./2.0;
			t24 = riw2./2.0;
			t25 = riw3./2.0;
			t26 = riw4./2.0;
			t27 = riw5./2.0;
			t28 = riw6./2.0;
			t29 = cos(t2);
			t30 = cos(t3);
			t31 = cos(t4);
			t32 = cos(t11);
			t33 = cos(t12);
			t34 = cos(t13);
			t35 = cos(t14);
			t36 = cos(t15);
			t37 = cos(t16);
			t38 = cos(t17);
			t39 = cos(t18);
			t40 = cos(t19);
			t41 = cos(t20);
			t42 = cos(t21);
			t43 = cos(t22);
			t44 = sin(t2);
			t45 = sin(t3);
			t46 = sin(t4);
			t47 = sin(t11);
			t48 = sin(t12);
			t49 = sin(t13);
			t50 = sin(t14);
			t51 = sin(t15);
			t52 = sin(t16);
			t53 = sin(t17);
			t54 = sin(t18);
			t55 = sin(t19);
			t56 = sin(t20);
			t57 = sin(t21);
			t58 = sin(t22);
			t59 = t5+t23;
			t60 = t6+t24;
			t61 = t7+t25;
			t62 = t8+t26;
			t63 = t9+t27;
			t64 = t10+t28;
			t65 = t34.*t36.*t52;
			t66 = t34.*t37.*t51;
			t67 = t36.*t37.*t49;
			t68 = t34.*t51.*t52;
			t69 = t36.*t49.*t52;
			t70 = t37.*t49.*t51;
			t71 = t49.*t51.*t52;
			t72 = t34.*t36.*t37;
			t73 = t71+t72;
			mt1 = [1.0,0.0,0.0,0.0,t32.*t33.*t35+t47.*t48.*t50,t33.*t35.*t47-t32.*t48.*t50,t32.*t35.*t48+t33.*t47.*t50,t32.*t33.*t50-t35.*t47.*t48,t41.*t42.*t43+t56.*t57.*t58,t42.*t43.*t56-t41.*t57.*t58,t41.*t43.*t57+t42.*t56.*t58,t41.*t42.*t58-t43.*t56.*t57,t29.*t30.*t31+t44.*t45.*t46,t29.*t31.*t45-t30.*t44.*t46,t29.*t30.*t46+t31.*t44.*t45,t30.*t31.*t44-t29.*t45.*t46,1.0,0.0,0.0,0.0,t73,t67-t68,t66+t69,t65-t70,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t38.*t39.*t40+t53.*t54.*t55,t39.*t40.*t53-t38.*t54.*t55,t38.*t40.*t54+t39.*t53.*t55];
			mt2 = [t38.*t39.*t55-t40.*t53.*t54,cos(t59),0.0,0.0,sin(t59),cos(t60),sin(t60),0.0,0.0,cos(t61),sin(t61),0.0,0.0,cos(t62),sin(t62),0.0,0.0,cos(t63),0.0,sin(t63),0.0,cos(t64),0.0,0.0,sin(t64),t73,-t67+t68,-t66-t69,-t65+t70,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0];
			out1 = reshape([mt1,mt2],4,26);
			end
			function out1 = linkCorAcc(t,in2,in3,in4,in5)
			%linkCorAcc
			%    OUT1 = linkCorAcc(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:19:06
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = sin(q1dev__dt_0_);
			t5 = sin(q1sup__dt_0_);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-q1dev__dt_1_.*q1sup__dt_1_.*t2-q1dev__dt_1_.*q1flex__dt_1_.*t3.*t4-q1flex__dt_1_.*q1sup__dt_1_.*t2.*t5,-q1dev__dt_1_.*q1sup__dt_1_.*t4+q1dev__dt_1_.*q1flex__dt_1_.*t2.*t3-q1flex__dt_1_.*q1sup__dt_1_.*t4.*t5,-q1flex__dt_1_.*q1sup__dt_1_.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,26]);
			end
			function out1 = linkJacobian(t,in2,in3,in4,in5)
			%linkJacobian
			%    OUT1 = linkJacobian(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:19:07
			q1dev__dt_0_ = in2(1,:);
			q1sup__dt_0_ = in2(3,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1sup__dt_0_);
			t4 = sin(q1dev__dt_0_);
			mt1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,t2.*t3,t3.*t4,-sin(q1sup__dt_0_),0.0,0.0,0.0,-t4,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt3 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt4 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0];
			mt5 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt6 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt7 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			mt8 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			out1 = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8],6,234);
			end
		end
		function [inertM,inertGF,inertCenMat,inertCorMatLeft,inertCorMatRight]=Inertial_TAWE(t,q,p,u,s,TF_Global,Vel_Global,Cor_Global,Jac_Global,Full_Flag_)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    BaseFrameList=[6   7  21  22  23  24  25  26];
		    MassList=Mass_TAWE(t,q,p,s,u);
		    if Full_Flag_
		        MomentList=reshape(Moment_TAWE(t,q,p,s,u),[3 3 numel(MassList)]);
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
			function out1 = Mass_TAWE(t,in2,in3,in4,in5)
			%Mass_TAWE
			%    OUT1 = Mass_TAWE(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:05
			m12 = in3(112,:);
			m13 = in3(63,:);
			m22 = in3(70,:);
			m23 = in3(77,:);
			m24 = in3(84,:);
			m25 = in3(91,:);
			m26 = in3(98,:);
			m27 = in3(105,:);
			out1 = [m12,m13,m22,m23,m24,m25,m26,m27];
			end
			function out1 = Moment_TAWE(t,in2,in3,in4,in5)
			%Moment_TAWE
			%    OUT1 = Moment_TAWE(T,IN2,IN3,IN4,IN5)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:05
			i12xx = in3(113,:);
			i12xy = in3(114,:);
			i13xx = in3(64,:);
			i12xz = in3(115,:);
			i12yy = in3(116,:);
			i13xy = in3(65,:);
			i12yz = in3(117,:);
			i13xz = in3(66,:);
			i13yy = in3(67,:);
			i12zz = in3(118,:);
			i13yz = in3(68,:);
			i13zz = in3(69,:);
			i22xx = in3(71,:);
			i22xy = in3(72,:);
			i23xx = in3(78,:);
			i22xz = in3(73,:);
			i22yy = in3(74,:);
			i23xy = in3(79,:);
			i24xx = in3(85,:);
			i22yz = in3(75,:);
			i23xz = in3(80,:);
			i23yy = in3(81,:);
			i24xy = in3(86,:);
			i25xx = in3(92,:);
			i22zz = in3(76,:);
			i23yz = in3(82,:);
			i24xz = in3(87,:);
			i24yy = in3(88,:);
			i25xy = in3(93,:);
			i26xx = in3(99,:);
			i23zz = in3(83,:);
			i24yz = in3(89,:);
			i25xz = in3(94,:);
			i25yy = in3(95,:);
			i26xy = in3(100,:);
			i27xx = in3(106,:);
			i24zz = in3(90,:);
			i25yz = in3(96,:);
			i26xz = in3(101,:);
			i26yy = in3(102,:);
			i27xy = in3(107,:);
			i25zz = in3(97,:);
			i26yz = in3(103,:);
			i27xz = in3(108,:);
			i27yy = in3(109,:);
			i26zz = in3(104,:);
			i27yz = in3(110,:);
			i27zz = in3(111,:);
			out1 = reshape([i12xx,i12xy,i12xz,i12xy,i12yy,i12yz,i12xz,i12yz,i12zz,i13xx,i13xy,i13xz,i13xy,i13yy,i13yz,i13xz,i13yz,i13zz,i22xx,i22xy,i22xz,i22xy,i22yy,i22yz,i22xz,i22yz,i22zz,i23xx,i23xy,i23xz,i23xy,i23yy,i23yz,i23xz,i23yz,i23zz,i24xx,i24xy,i24xz,i24xy,i24yy,i24yz,i24xz,i24yz,i24zz,i25xx,i25xy,i25xz,i25xy,i25yy,i25yz,i25xz,i25yz,i25zz,i26xx,i26xy,i26xz,i26xy,i26yy,i26yz,i26xz,i26yz,i26zz,i27xx,i27xy,i27xz,i27xy,i27yy,i27yz,i27xz,i27yz,i27zz],[3,24]);
			end
		end
		function [CollectGF]=Force_TAWE(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		 
		    CollectGF=zeros(numel(q)/2,1);
		    ForceNum=9;
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
				        SubFrame=[0];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ForceFrameVecJac=ForceFVJac_2(t,q,p,u,s,SubSubs);
				        SubJac=ForceJac_2(t,q,p,u,s,SubSubs);
				        SubEff=ForceEff_2(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubJac=SubJac+ForceFrameVecJac(:,:,sfCount)*JacSubs(:,:,sfCount);
				            end
				        end
				        SubGF=SubJac.'*SubEff;
				    case 3
				        SubFrame=[7];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[7];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_3(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				    case 4
				        SubFrame=[21];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[21];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_4(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				    case 5
				        SubFrame=[22];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[22];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_5(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				    case 6
				        SubFrame=[23];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[23];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_6(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				    case 7
				        SubFrame=[24];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[24];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_7(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				    case 8
				        SubFrame=[25];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[25];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_8(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				    case 9
				        SubFrame=[26];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[26];
				        RefFrame=[1];
				        SubJac=Jac_Global(1:3,:,ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*ForceEff_9(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				end
		        CollectGF(:,FCount)=SubGF;
		    end
			function out1 = ForceJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceJac_1
			%    OUT1 = ForceJac_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:05
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,9]);
			end
			function out1 = ForceJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceJac_2
			%    OUT1 = ForceJac_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:05
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0],[6,9]);
			end
			function out1 = ForceJac_3(t,in2,in3,in4,in5,in6)
			%ForceJac_3
			%    OUT1 = ForceJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			out1 = 0.0;
			end
			function out1 = ForceJac_4(t,in2,in3,in4,in5,in6)
			%ForceJac_4
			%    OUT1 = ForceJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			out1 = 0.0;
			end
			function out1 = ForceJac_5(t,in2,in3,in4,in5,in6)
			%ForceJac_5
			%    OUT1 = ForceJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			out1 = 0.0;
			end
			function out1 = ForceJac_6(t,in2,in3,in4,in5,in6)
			%ForceJac_6
			%    OUT1 = ForceJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			out1 = 0.0;
			end
			function out1 = ForceJac_7(t,in2,in3,in4,in5,in6)
			%ForceJac_7
			%    OUT1 = ForceJac_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			out1 = 0.0;
			end
			function out1 = ForceJac_8(t,in2,in3,in4,in5,in6)
			%ForceJac_8
			%    OUT1 = ForceJac_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			out1 = 0.0;
			end
			function out1 = ForceJac_9(t,in2,in3,in4,in5,in6)
			%ForceJac_9
			%    OUT1 = ForceJac_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			out1 = 0.0;
			end
			function out1 = ForceEff_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceEff_1
			%    OUT1 = ForceEff_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:05
			b1act = in3(4,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_1_ = in2(12,:);
			out1 = [-b1act.*q1dev__dt_1_;-b1act.*q1flex__dt_1_;-b1act.*q1sup__dt_1_];
			end
			function out1 = ForceEff_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceEff_2
			%    OUT1 = ForceEff_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			b2act = in3(5,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_1_ = in2(18,:);
			out1 = [-b2act.*q2w1__dt_1_;-b2act.*q2w2__dt_1_;-b2act.*q2w3__dt_1_;-b2act.*q2w4__dt_1_;-b2act.*q2w5__dt_1_;-b2act.*q2w6__dt_1_];
			end
			function out1 = ForceEff_3(t,in2,in3,in4,in5,in6)
			%ForceEff_3
			%    OUT1 = ForceEff_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			gUncertain = in3(2,:);
			m13 = in3(63,:);
			out1 = [0.0;0.0;-gUncertain.*m13];
			end
			function out1 = ForceEff_4(t,in2,in3,in4,in5,in6)
			%ForceEff_4
			%    OUT1 = ForceEff_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			g = in3(1,:);
			m22 = in3(70,:);
			out1 = [0.0;0.0;-g.*m22];
			end
			function out1 = ForceEff_5(t,in2,in3,in4,in5,in6)
			%ForceEff_5
			%    OUT1 = ForceEff_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			g = in3(1,:);
			m23 = in3(77,:);
			out1 = [0.0;0.0;-g.*m23];
			end
			function out1 = ForceEff_6(t,in2,in3,in4,in5,in6)
			%ForceEff_6
			%    OUT1 = ForceEff_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			g = in3(1,:);
			m24 = in3(84,:);
			out1 = [0.0;0.0;-g.*m24];
			end
			function out1 = ForceEff_7(t,in2,in3,in4,in5,in6)
			%ForceEff_7
			%    OUT1 = ForceEff_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			g = in3(1,:);
			m25 = in3(91,:);
			out1 = [0.0;0.0;-g.*m25];
			end
			function out1 = ForceEff_8(t,in2,in3,in4,in5,in6)
			%ForceEff_8
			%    OUT1 = ForceEff_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			g = in3(1,:);
			m26 = in3(98,:);
			out1 = [0.0;0.0;-g.*m26];
			end
			function out1 = ForceEff_9(t,in2,in3,in4,in5,in6)
			%ForceEff_9
			%    OUT1 = ForceEff_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			g = in3(1,:);
			m27 = in3(105,:);
			out1 = [0.0;0.0;-g.*m27];
			end
			function out1 = ForceFVJac_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceFVJac_1
			%    OUT1 = ForceFVJac_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:05
			out1 = 0.0;
			end
			function out1 = ForceFVJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ForceFVJac_2
			%    OUT1 = ForceFVJac_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			out1 = 0.0;
			end
			function out1 = ForceFVJac_3(t,in2,in3,in4,in5,in6)
			%ForceFVJac_3
			%    OUT1 = ForceFVJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = ForceFVJac_4(t,in2,in3,in4,in5,in6)
			%ForceFVJac_4
			%    OUT1 = ForceFVJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:06
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = ForceFVJac_5(t,in2,in3,in4,in5,in6)
			%ForceFVJac_5
			%    OUT1 = ForceFVJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = ForceFVJac_6(t,in2,in3,in4,in5,in6)
			%ForceFVJac_6
			%    OUT1 = ForceFVJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = ForceFVJac_7(t,in2,in3,in4,in5,in6)
			%ForceFVJac_7
			%    OUT1 = ForceFVJac_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:07
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = ForceFVJac_8(t,in2,in3,in4,in5,in6)
			%ForceFVJac_8
			%    OUT1 = ForceFVJac_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = ForceFVJac_9(t,in2,in3,in4,in5,in6)
			%ForceFVJac_9
			%    OUT1 = ForceFVJac_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
		end
		function [CollectGF,CollectInputJac]=Input_TAWE(t,q,p,u,s,TF_Global,TransDis_Global,Vel_Global,Jac_Global,Quat_Global)
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    CollectGF=zeros(numel(q)/2,1);
		    CollectInputJac=zeros(numel(q)/2,numel(u));
		    ForceNum=9;
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
				    case 7
				        SubFrame=[21];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[-21];
				        RefFrame=[14];
				        SubJac=Jac_Global(4:6,:,-ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_7(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_7(t,q,p,u,s,SubSubs);
				    case 8
				        SubFrame=[21];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[-21];
				        RefFrame=[15];
				        SubJac=Jac_Global(4:6,:,-ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_8(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_8(t,q,p,u,s,SubSubs);
				    case 9
				        SubFrame=[22];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ActFrame=[-22];
				        RefFrame=[15];
				        SubJac=Jac_Global(4:6,:,-ActFrame);
				        SubEff=TF_Global(1:3,1:3,RefFrame)*InputEff_9(t,q,p,u,s,SubSubs);
				        SubGF=SubJac.'*SubEff;
				        SubInputJac=SubJac.'*InputEffJac_9(t,q,p,u,s,SubSubs);
				end
		        CollectGF(:,FCount)=SubGF;
		        CollectInputJac(:,:,FCount)=SubInputJac;
		    end
			function out1 = ForceJac_1(t,in2,in3,in4,in5,in6)
			%ForceJac_1
			%    OUT1 = ForceJac_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			out1 = 0.0;
			end
			function out1 = ForceJac_2(t,in2,in3,in4,in5,in6)
			%ForceJac_2
			%    OUT1 = ForceJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			out1 = 0.0;
			end
			function out1 = ForceJac_3(t,in2,in3,in4,in5,in6)
			%ForceJac_3
			%    OUT1 = ForceJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			out1 = 0.0;
			end
			function out1 = ForceJac_4(t,in2,in3,in4,in5,in6)
			%ForceJac_4
			%    OUT1 = ForceJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			out1 = 0.0;
			end
			function out1 = ForceJac_5(t,in2,in3,in4,in5,in6)
			%ForceJac_5
			%    OUT1 = ForceJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			out1 = 0.0;
			end
			function out1 = ForceJac_6(t,in2,in3,in4,in5,in6)
			%ForceJac_6
			%    OUT1 = ForceJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:11
			out1 = 0.0;
			end
			function out1 = ForceJac_7(t,in2,in3,in4,in5,in6)
			%ForceJac_7
			%    OUT1 = ForceJac_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:11
			out1 = 0.0;
			end
			function out1 = ForceJac_8(t,in2,in3,in4,in5,in6)
			%ForceJac_8
			%    OUT1 = ForceJac_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:12
			out1 = 0.0;
			end
			function out1 = ForceJac_9(t,in2,in3,in4,in5,in6)
			%ForceJac_9
			%    OUT1 = ForceJac_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:13
			out1 = 0.0;
			end
			function out1 = InputEff_1(t,in2,in3,in4,in5,in6)
			%InputEff_1
			%    OUT1 = InputEff_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:08
			g = in3(1,:);
			handmpt1 = in4(5,:);
			out1 = [0.0;0.0;-g.*handmpt1];
			end
			function out1 = InputEff_2(t,in2,in3,in4,in5,in6)
			%InputEff_2
			%    OUT1 = InputEff_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			g = in3(1,:);
			handmpt2 = in4(6,:);
			out1 = [0.0;0.0;-g.*handmpt2];
			end
			function out1 = InputEff_3(t,in2,in3,in4,in5,in6)
			%InputEff_3
			%    OUT1 = InputEff_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			g = in3(1,:);
			handmpt3 = in4(7,:);
			out1 = [0.0;0.0;-g.*handmpt3];
			end
			function out1 = InputEff_4(t,in2,in3,in4,in5,in6)
			%InputEff_4
			%    OUT1 = InputEff_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			g = in3(1,:);
			handmpt4 = in4(8,:);
			out1 = [0.0;0.0;-g.*handmpt4];
			end
			function out1 = InputEff_5(t,in2,in3,in4,in5,in6)
			%InputEff_5
			%    OUT1 = InputEff_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			u1dev = in4(1,:);
			out1 = [0.0;0.0;u1dev];
			end
			function out1 = InputEff_6(t,in2,in3,in4,in5,in6)
			%InputEff_6
			%    OUT1 = InputEff_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:11
			u1flex = in4(2,:);
			out1 = [u1flex;0.0;0.0];
			end
			function out1 = InputEff_7(t,in2,in3,in4,in5,in6)
			%InputEff_7
			%    OUT1 = InputEff_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:12
			u21 = in4(3,:);
			out1 = [0.0;0.0;u21];
			end
			function out1 = InputEff_8(t,in2,in3,in4,in5,in6)
			%InputEff_8
			%    OUT1 = InputEff_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:12
			u22 = in4(4,:);
			out1 = [-u22;0.0;0.0];
			end
			function out1 = InputEff_9(t,in2,in3,in4,in5,in6)
			%InputEff_9
			%    OUT1 = InputEff_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:13
			u22 = in4(4,:);
			out1 = [u22;0.0;0.0];
			end
			function out1 = InputEffJac_1(t,in2,in3,in4,in5,in6)
			%InputEffJac_1
			%    OUT1 = InputEffJac_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_2(t,in2,in3,in4,in5,in6)
			%InputEffJac_2
			%    OUT1 = InputEffJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_3(t,in2,in3,in4,in5,in6)
			%InputEffJac_3
			%    OUT1 = InputEffJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_4(t,in2,in3,in4,in5,in6)
			%InputEffJac_4
			%    OUT1 = InputEffJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			g = in3(1,:);
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g],[3,8]);
			end
			function out1 = InputEffJac_5(t,in2,in3,in4,in5,in6)
			%InputEffJac_5
			%    OUT1 = InputEffJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU1z = in3(26,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
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
			out1 = reshape([t27.*(t2.*t10-t4.*t8.*t9)+t5.*t6.*(t8.*t10+t2.*t4.*t9)+t3.*t4.*t26,-t27.*(t23.*-2.0-t24.*2.0+t23.*t24.*4.0+t19.*t20.*cos(t17).*sin(t16).*sin(t17).*sin(t18).*8.0+1.0)-t5.*t6.*(t4.*t8-t2.*t9.*t10)+t3.*t10.*t26,-t9.*t26-t3.*t8.*t27+t2.*t3.*t5.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_6(t,in2,in3,in4,in5,in6)
			%InputEffJac_6
			%    OUT1 = InputEffJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:11
			q1dev__dt_0_ = in2(1,:);
			q1sup__dt_0_ = in2(3,:);
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU1z = in3(26,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
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
			out1 = reshape([0.0,0.0,0.0,-t62.*(t12.*t14+t4.*t6.*t13)-t67.*(t4.*t14-t6.*t12.*t13)-t5.*t6.*t64,t62.*(t6.*t12-t4.*t13.*t14)+t67.*(t40.*-2.0-t41.*2.0+t40.*t41.*4.0+t27.*t28.*cos(t22).*sin(t21).*sin(t22).*sin(t23).*8.0+1.0)-t5.*t14.*t64,t13.*t64-t4.*t5.*t62+t5.*t12.*t67,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_7(t,in2,in3,in4,in5,in6)
			%InputEffJac_7
			%    OUT1 = InputEffJac_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:12
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU1z = in3(26,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			t2 = cos(rf1IMU1x);
			t3 = cos(rf1IMU1y);
			t4 = cos(rf1IMU1z);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = sin(rf1IMU1x);
			t9 = sin(rf1IMU1y);
			t10 = sin(rf1IMU1z);
			t11 = sin(rf2Medx);
			t12 = sin(rf2Medy);
			t13 = sin(rf2Medz);
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
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t27.*(t2.*t10-t4.*t8.*t9)+t5.*t6.*(t8.*t10+t2.*t4.*t9)+t3.*t4.*t26,-t27.*(t23.*-2.0-t24.*2.0+t23.*t24.*4.0+t19.*t20.*cos(t17).*sin(t16).*sin(t17).*sin(t18).*8.0+1.0)-t5.*t6.*(t4.*t8-t2.*t9.*t10)+t3.*t10.*t26,-t9.*t26-t3.*t8.*t27+t2.*t3.*t5.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_8(t,in2,in3,in4,in5,in6)
			%InputEffJac_8
			%    OUT1 = InputEffJac_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:12
			q2w1__dt_0_ = in2(4,:);
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU1z = in3(26,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			riw1 = in3(36,:);
			t2 = cos(rf1IMU1x);
			t3 = cos(rf1IMU1y);
			t4 = cos(rf1IMU1z);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = sin(rf1IMU1x);
			t9 = sin(rf1IMU1y);
			t10 = sin(rf1IMU1z);
			t11 = sin(rf2Medx);
			t12 = sin(rf2Medy);
			t13 = sin(rf2Medz);
			t14 = q2w1__dt_0_+riw1;
			t18 = rf1IMU1x./2.0;
			t19 = rf1IMU1y./2.0;
			t20 = rf1IMU1z./2.0;
			t21 = rf2Medx./2.0;
			t22 = rf2Medy./2.0;
			t23 = rf2Medz./2.0;
			t15 = cos(t14);
			t16 = sin(t14);
			t17 = t5.*t13;
			t24 = cos(t18);
			t25 = cos(t20);
			t26 = cos(t21);
			t27 = cos(t22);
			t28 = cos(t23);
			t29 = sin(t21);
			t30 = sin(t22);
			t31 = sin(t23);
			t32 = t7.*t11.*t12;
			t33 = t12.*t15;
			t34 = t24.^2;
			t35 = t25.^2;
			t36 = t26.^2;
			t37 = t28.^2;
			t38 = -t32;
			t39 = t6.*t7.*t15;
			t40 = t6.*t13.*t15;
			t41 = t6.*t11.*t16;
			t53 = t26.*t27.*t28.*t29.*t30.*t31.*8.0;
			t42 = t36.*2.0;
			t43 = t37.*2.0;
			t44 = -t39;
			t45 = -t41;
			t48 = t17+t38;
			t49 = t36.*t37.*4.0;
			t46 = -t42;
			t47 = -t43;
			t50 = t33+t45;
			t51 = t16.*t48;
			t52 = t44+t51;
			t54 = t46+t47+t49+t53+1.0;
			t55 = t16.*t54;
			t56 = t40+t55;
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t50.*(t8.*t10+t2.*t4.*t9)+t56.*(t2.*t10-t4.*t8.*t9)-t3.*t4.*(t39-t51),-t50.*(t4.*t8-t2.*t9.*t10)-t56.*(t34.*-2.0-t35.*2.0+t34.*t35.*4.0+t24.*t25.*cos(t19).*sin(t18).*sin(t19).*sin(t20).*8.0+1.0)-t3.*t10.*(t39-t51),t9.*(t39-t51)+t2.*t3.*t50-t3.*t8.*t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputEffJac_9(t,in2,in3,in4,in5,in6)
			%InputEffJac_9
			%    OUT1 = InputEffJac_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:13
			q2w1__dt_0_ = in2(4,:);
			rf1IMU1x = in3(24,:);
			rf1IMU1y = in3(25,:);
			rf1IMU1z = in3(26,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			riw1 = in3(36,:);
			t2 = cos(rf1IMU1x);
			t3 = cos(rf1IMU1y);
			t4 = cos(rf1IMU1z);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = sin(rf1IMU1x);
			t9 = sin(rf1IMU1y);
			t10 = sin(rf1IMU1z);
			t11 = sin(rf2Medx);
			t12 = sin(rf2Medy);
			t13 = sin(rf2Medz);
			t14 = q2w1__dt_0_+riw1;
			t18 = rf1IMU1x./2.0;
			t19 = rf1IMU1y./2.0;
			t20 = rf1IMU1z./2.0;
			t21 = rf2Medx./2.0;
			t22 = rf2Medy./2.0;
			t23 = rf2Medz./2.0;
			t15 = cos(t14);
			t16 = sin(t14);
			t17 = t5.*t13;
			t24 = cos(t18);
			t25 = cos(t20);
			t26 = cos(t21);
			t27 = cos(t22);
			t28 = cos(t23);
			t29 = sin(t21);
			t30 = sin(t22);
			t31 = sin(t23);
			t32 = t7.*t11.*t12;
			t33 = t12.*t15;
			t34 = t24.^2;
			t35 = t25.^2;
			t36 = t26.^2;
			t37 = t28.^2;
			t38 = -t32;
			t39 = t6.*t7.*t15;
			t40 = t6.*t13.*t15;
			t41 = t6.*t11.*t16;
			t53 = t26.*t27.*t28.*t29.*t30.*t31.*8.0;
			t42 = t36.*2.0;
			t43 = t37.*2.0;
			t44 = -t39;
			t45 = -t41;
			t48 = t17+t38;
			t49 = t36.*t37.*4.0;
			t46 = -t42;
			t47 = -t43;
			t50 = t33+t45;
			t51 = t16.*t48;
			t52 = t44+t51;
			t54 = t46+t47+t49+t53+1.0;
			t55 = t16.*t54;
			t56 = t40+t55;
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t50.*(t8.*t10+t2.*t4.*t9)-t56.*(t2.*t10-t4.*t8.*t9)+t3.*t4.*(t39-t51),t50.*(t4.*t8-t2.*t9.*t10)+t56.*(t34.*-2.0-t35.*2.0+t34.*t35.*4.0+t24.*t25.*cos(t19).*sin(t18).*sin(t19).*sin(t20).*8.0+1.0)+t3.*t10.*(t39-t51),-t9.*(t39-t51)-t2.*t3.*t50+t3.*t8.*t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,8]);
			end
			function out1 = InputFVJac_1(t,in2,in3,in4,in5,in6)
			%InputFVJac_1
			%    OUT1 = InputFVJac_1(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_2(t,in2,in3,in4,in5,in6)
			%InputFVJac_2
			%    OUT1 = InputFVJac_2(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:09
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_3(t,in2,in3,in4,in5,in6)
			%InputFVJac_3
			%    OUT1 = InputFVJac_3(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_4(t,in2,in3,in4,in5,in6)
			%InputFVJac_4
			%    OUT1 = InputFVJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:10
			out1 = reshape([1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
			end
			function out1 = InputFVJac_5(t,in2,in3,in4,in5,in6)
			%InputFVJac_5
			%    OUT1 = InputFVJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:11
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
			function out1 = InputFVJac_6(t,in2,in3,in4,in5,in6)
			%InputFVJac_6
			%    OUT1 = InputFVJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:11
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
			function out1 = InputFVJac_7(t,in2,in3,in4,in5,in6)
			%InputFVJac_7
			%    OUT1 = InputFVJac_7(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:12
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
			function out1 = InputFVJac_8(t,in2,in3,in4,in5,in6)
			%InputFVJac_8
			%    OUT1 = InputFVJac_8(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:13
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
			function out1 = InputFVJac_9(t,in2,in3,in4,in5,in6)
			%InputFVJac_9
			%    OUT1 = InputFVJac_9(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:13
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0],[3,6]);
			end
		end
		function [ConsJac,ConsCor,ConsGF]=Constraint_TAWE(t,q,p,u,s,TransDis_Global,Vel_Global,Cor_Global,Jac_Global,Quat_Global)
		%% Constraint dynamic property calculator (Toolbox Internal Use)
		%% 
		%% Author: Jiamin Wang; Updated: 2021-12-15;
		    ConsNum=7;
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
				    case 2
				        SubFrame=[0];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_2(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_2(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_2(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_2(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            CorSubs=Cor_Global(:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
				                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
				            end
				        end
				    case 3
				        SubFrame=[0];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_3(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_3(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_3(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_3(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            CorSubs=Cor_Global(:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
				                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
				            end
				        end
				    case 4
				        SubFrame=[5  20];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_4(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_4(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_4(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_4(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            CorSubs=Cor_Global(:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
				                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
				            end
				        end
				    case 5
				        SubFrame=[5  20];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_5(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_5(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_5(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_5(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            CorSubs=Cor_Global(:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
				                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
				            end
				        end
				    case 6
				        SubFrame=[5  20];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_6(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_6(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_6(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_6(t,q,p,u,s,SubSubs);
				        if SubFrame(1)~=0
				            JacSubs=Jac_Global(:,:,SubFrame);
				            CorSubs=Cor_Global(:,SubFrame);
				            for sfCount=1:numel(SubFrame)
				                SubConsJac=SubConsJac+ConsFrameVecJac(sfCount,:)*JacSubs(:,:,sfCount);
				                SubConsCor=SubConsCor+ConsFrameVecJac(sfCount,:)*CorSubs(:,sfCount);
				            end
				        end
				    case 7
				        SubFrame=[0];
				        SubSubs=0;
				        if SubFrame(1)~=0
				            SubSubs=[TransDis_Global(:,SubFrame);Vel_Global(:,SubFrame);Quat_Global(:,SubFrame)];
				        end
				        ConsFrameVecJac=ConsFVJac_7(t,q,p,u,s,SubSubs);
				        SubConsJac=ConsJac_7(t,q,p,u,s,SubSubs);
				        SubConsCor=ConsCor_7(t,q,p,u,s,SubSubs);
				        SubConsGF=ConsGF_7(t,q,p,u,s,SubSubs);
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
			%    14-Apr-2023 16:21:14
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			q2w1__dt_0_ = in2(4,:);
			q2w2__dt_0_ = in2(5,:);
			q2w3__dt_0_ = in2(6,:);
			q2w4__dt_0_ = in2(7,:);
			q2w5__dt_0_ = in2(8,:);
			q2w6__dt_0_ = in2(9,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = cos(rf1Devx);
			t9 = cos(rf1Devy);
			t10 = cos(rf1Devz);
			t11 = sin(q1dev__dt_0_);
			t12 = sin(q1flex__dt_0_);
			t13 = sin(q1sup__dt_0_);
			t14 = sin(rf2Medx);
			t15 = sin(rf2Medy);
			t16 = sin(rf2Medz);
			t17 = sin(rf1Devx);
			t18 = sin(rf1Devy);
			t19 = sin(rf1Devz);
			t20 = l25+l13z;
			t21 = q2w1__dt_0_+riw1;
			t22 = q2w2__dt_0_+riw2;
			t23 = q2w3__dt_0_+riw3;
			t24 = q2w4__dt_0_+riw4;
			t25 = q2w5__dt_0_+riw5;
			t26 = q2w6__dt_0_+riw6;
			t45 = q1dev__dt_0_./2.0;
			t46 = q1flex__dt_0_./2.0;
			t47 = q1sup__dt_0_./2.0;
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = cos(t23);
			t30 = cos(t24);
			t31 = cos(t25);
			t32 = cos(t26);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = sin(t23);
			t36 = sin(t24);
			t37 = sin(t25);
			t38 = sin(t26);
			t39 = t2.*t3;
			t40 = t5.*t16;
			t41 = t8.*t19;
			t42 = t11.*t12;
			t43 = t14.*t16;
			t44 = t17.*t19;
			t48 = cos(t45);
			t49 = cos(t46);
			t50 = cos(t47);
			t51 = sin(t45);
			t52 = sin(t46);
			t53 = sin(t47);
			t58 = t7.*t14.*t15;
			t59 = t10.*t17.*t18;
			t62 = t5.*t7.*t15;
			t63 = t8.*t10.*t18;
			t54 = l13x.*t32;
			t55 = l13y.*t32;
			t56 = l13x.*t38;
			t57 = l13y.*t38;
			t60 = t13.*t42;
			t61 = t13.*t39;
			t64 = t20.*t31;
			t65 = t20.*t37;
			t67 = t48.^2;
			t68 = t49.^2;
			t69 = -t58;
			t70 = -t59;
			t74 = t43+t62;
			t75 = t44+t63;
			t66 = -t57;
			t71 = t55+t56;
			t72 = t39+t60;
			t73 = t42+t61;
			t78 = t40+t69;
			t79 = t41+t70;
			t76 = l23+t71;
			t77 = t54+t66;
			t88 = t30.*t37.*t71;
			t89 = t36.*t37.*t71;
			t80 = l24+t77;
			t81 = t30.*t76;
			t82 = t36.*t76;
			t83 = t30.*t77;
			t84 = t36.*t77;
			t91 = -t89;
			t85 = -t81;
			t86 = t31.*t80;
			t87 = t37.*t80;
			t102 = t84+t88;
			t103 = t83+t91;
			t90 = -t87;
			t92 = t65+t86;
			t104 = t29.*t102;
			t105 = t35.*t102;
			t106 = t29.*t103;
			t107 = t35.*t103;
			t93 = t64+t90;
			t96 = t29.*t30.*t92;
			t97 = t29.*t36.*t92;
			t98 = t30.*t35.*t92;
			t99 = t35.*t36.*t92;
			t108 = -t105;
			t124 = t104+t107;
			t127 = -t28.*(t105-t106);
			t94 = t30.*t93;
			t95 = t36.*t93;
			t101 = -t99;
			t120 = t97+t98;
			t125 = t106+t108;
			t126 = t34.*t124;
			t100 = -t95;
			t109 = t82+t94;
			t110 = t85+t95;
			t121 = t96+t101;
			t122 = t28.*t120;
			t128 = -t126;
			t111 = l22+t81+t100;
			t112 = t29.*t109;
			t113 = t35.*t109;
			t114 = -t29.*(t81+t100);
			t115 = -t35.*(t81+t100);
			t119 = t35.*(t81+t100);
			t123 = t34.*t121;
			t144 = t127+t128;
			t116 = -t113;
			t117 = t29.*t111;
			t118 = t35.*t111;
			t129 = t113+t114;
			t131 = t112+t119;
			t143 = t122+t123;
			t130 = t112+t118;
			t132 = t116+t117;
			t134 = t34.*t129;
			t137 = t28.*t131;
			t139 = -t34.*(t113-t117);
			t133 = l21+t132;
			t135 = t28.*t130;
			t136 = t34.*t130;
			t138 = -t134;
			t140 = t28.*t133;
			t141 = t34.*t133;
			t145 = t137+t138;
			t146 = t135+t139;
			t142 = -t140;
			t147 = t135+t141;
			t148 = t136+t142;
			mt1 = [-t79.*(l12z.*t73+l12y.*(t48.*t51.*2.0-t48.*t51.*t68.*4.0-t49.*t50.*t51.^2.*t52.*t53.*4.0+t49.*t50.*t52.*t53.*t67.*4.0)+l12x.*t2.*t4)-t9.*t10.*(-l12z.*(t2.*t12-t3.*t11.*t13)+l12y.*t72+l12x.*t4.*t11),t79.*(l12z.*t72-l12y.*(t49.*t52.*2.0-t49.*t52.*t67.*4.0-t48.*t50.*t51.*t52.^2.*t53.*4.0+t48.*t50.*t51.*t53.*t68.*4.0))+t75.*(l12y.*t3.*t4-l12z.*t4.*t12)+t9.*t10.*(l12z.*(t3.*t11-t2.*t12.*t13)+l12y.*t73)];
			mt2 = [-t79.*(l12y.*(t48.*t49.*t50.^2.*t51.*t52.*4.0-t48.*t49.*t51.*t52.*t53.^2.*4.0)-l12x.*t11.*t13+l12z.*t3.*t4.*t11)-t75.*(l12x.*t4+l12z.*t3.*t13+l12y.*t12.*t13)+t9.*t10.*(-l12x.*t2.*t13+l12z.*t4.*t39+l12y.*t2.*t4.*t12),t78.*(t27.*t92+t33.*t148)+t6.*t7.*(t33.*t92-t27.*t148),t74.*t148-t27.*t78.*t147-t6.*t7.*t33.*t147,t74.*(t136+t28.*(t113-t117))-t27.*t78.*t146-t6.*t7.*t33.*t146,t74.*(t28.*t129+t34.*t131)+t27.*t78.*(t134-t137)+t6.*t7.*t33.*(t134-t137)];
			mt3 = [t78.*(t33.*t93+t27.*t143)+t74.*(t28.*t121-t34.*t120)-t6.*t7.*(t27.*t93-t33.*t143),-t78.*(t27.*(t126+t28.*(t105-t106))+t31.*t33.*t71)-t74.*(t28.*t124-t34.*(t105-t106))-t6.*t7.*(t33.*(t126+t28.*(t105-t106))-t27.*t31.*t71)];
			out1 = [mt1,mt2,mt3];
			end
			function out1 = ConsJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsJac_2
			%    OUT1 = ConsJac_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:18
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			q2w1__dt_0_ = in2(4,:);
			q2w2__dt_0_ = in2(5,:);
			q2w3__dt_0_ = in2(6,:);
			q2w4__dt_0_ = in2(7,:);
			q2w5__dt_0_ = in2(8,:);
			q2w6__dt_0_ = in2(9,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = cos(rf1Devx);
			t9 = cos(rf1Devy);
			t10 = cos(rf1Devz);
			t11 = sin(q1dev__dt_0_);
			t12 = sin(q1flex__dt_0_);
			t13 = sin(q1sup__dt_0_);
			t14 = sin(rf2Medx);
			t15 = sin(rf2Medy);
			t16 = sin(rf2Medz);
			t17 = sin(rf1Devx);
			t18 = sin(rf1Devy);
			t19 = sin(rf1Devz);
			t20 = l25+l13z;
			t21 = q2w1__dt_0_+riw1;
			t22 = q2w2__dt_0_+riw2;
			t23 = q2w3__dt_0_+riw3;
			t24 = q2w4__dt_0_+riw4;
			t25 = q2w5__dt_0_+riw5;
			t26 = q2w6__dt_0_+riw6;
			t43 = q1dev__dt_0_./2.0;
			t44 = q1flex__dt_0_./2.0;
			t45 = q1sup__dt_0_./2.0;
			t46 = rf2Medx./2.0;
			t47 = rf2Medy./2.0;
			t48 = rf2Medz./2.0;
			t49 = rf1Devx./2.0;
			t50 = rf1Devy./2.0;
			t51 = rf1Devz./2.0;
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = cos(t23);
			t30 = cos(t24);
			t31 = cos(t25);
			t32 = cos(t26);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = sin(t23);
			t36 = sin(t24);
			t37 = sin(t25);
			t38 = sin(t26);
			t39 = t2.*t3;
			t40 = t7.*t14;
			t41 = t10.*t17;
			t42 = t11.*t12;
			t52 = cos(t43);
			t53 = cos(t44);
			t54 = cos(t45);
			t55 = cos(t46);
			t56 = cos(t47);
			t57 = cos(t48);
			t58 = cos(t49);
			t59 = cos(t50);
			t60 = cos(t51);
			t61 = sin(t43);
			t62 = sin(t44);
			t63 = sin(t45);
			t64 = sin(t46);
			t65 = sin(t47);
			t66 = sin(t48);
			t67 = sin(t49);
			t68 = sin(t50);
			t69 = sin(t51);
			t74 = t5.*t15.*t16;
			t75 = t8.*t18.*t19;
			t70 = l13x.*t32;
			t71 = l13y.*t32;
			t72 = l13x.*t38;
			t73 = l13y.*t38;
			t76 = t13.*t42;
			t77 = t13.*t39;
			t78 = t20.*t31;
			t79 = t20.*t37;
			t81 = t52.^2;
			t82 = t53.^2;
			t83 = t55.^2;
			t84 = t57.^2;
			t85 = t58.^2;
			t86 = t60.^2;
			t87 = -t74;
			t88 = -t75;
			t127 = t55.*t56.*t57.*t64.*t65.*t66.*8.0;
			t128 = t58.*t59.*t60.*t67.*t68.*t69.*8.0;
			t80 = -t73;
			t89 = t83.*2.0;
			t90 = t84.*2.0;
			t91 = t85.*2.0;
			t92 = t86.*2.0;
			t97 = t71+t72;
			t98 = t39+t76;
			t99 = t42+t77;
			t102 = t40+t87;
			t103 = t41+t88;
			t105 = t83.*t84.*4.0;
			t106 = t85.*t86.*4.0;
			t93 = -t89;
			t94 = -t90;
			t95 = -t91;
			t96 = -t92;
			t100 = l23+t97;
			t101 = t70+t80;
			t114 = t30.*t37.*t97;
			t115 = t36.*t37.*t97;
			t104 = l24+t101;
			t107 = t30.*t100;
			t108 = t36.*t100;
			t109 = t30.*t101;
			t110 = t36.*t101;
			t117 = -t115;
			t149 = t93+t94+t105+t127+1.0;
			t150 = t95+t96+t106+t128+1.0;
			t111 = -t107;
			t112 = t31.*t104;
			t113 = t37.*t104;
			t130 = t110+t114;
			t131 = t109+t117;
			t116 = -t113;
			t118 = t79+t112;
			t132 = t29.*t130;
			t133 = t35.*t130;
			t134 = t29.*t131;
			t135 = t35.*t131;
			t119 = t78+t116;
			t122 = t29.*t30.*t118;
			t123 = t29.*t36.*t118;
			t124 = t30.*t35.*t118;
			t125 = t35.*t36.*t118;
			t136 = -t133;
			t154 = t132+t135;
			t157 = -t28.*(t133-t134);
			t120 = t30.*t119;
			t121 = t36.*t119;
			t129 = -t125;
			t148 = t123+t124;
			t155 = t134+t136;
			t156 = t34.*t154;
			t126 = -t121;
			t137 = t108+t120;
			t138 = t111+t121;
			t151 = t122+t129;
			t152 = t28.*t148;
			t158 = -t156;
			t139 = l22+t107+t126;
			t140 = t29.*t137;
			t141 = t35.*t137;
			t142 = -t29.*(t107+t126);
			t143 = -t35.*(t107+t126);
			t147 = t35.*(t107+t126);
			t153 = t34.*t151;
			t174 = t157+t158;
			t144 = -t141;
			t145 = t29.*t139;
			t146 = t35.*t139;
			t159 = t141+t142;
			t161 = t140+t147;
			t173 = t152+t153;
			t160 = t140+t146;
			t162 = t144+t145;
			t164 = t34.*t159;
			t167 = t28.*t161;
			t169 = -t34.*(t141-t145);
			t163 = l21+t162;
			t165 = t28.*t160;
			t166 = t34.*t160;
			t168 = -t164;
			t170 = t28.*t163;
			t171 = t34.*t163;
			t175 = t167+t168;
			t176 = t165+t169;
			t172 = -t170;
			t177 = t165+t171;
			t178 = t166+t172;
			mt1 = [t150.*(l12z.*t99+l12y.*(t52.*t61.*2.0-t52.*t61.*t82.*4.0-t53.*t54.*t61.^2.*t62.*t63.*4.0+t53.*t54.*t62.*t63.*t81.*4.0)+l12x.*t2.*t4)-t9.*t19.*(-l12z.*(t2.*t12-t3.*t11.*t13)+l12y.*t98+l12x.*t4.*t11),-t150.*(l12z.*t98-l12y.*(t53.*t62.*2.0-t53.*t62.*t81.*4.0-t52.*t54.*t61.*t62.^2.*t63.*4.0+t52.*t54.*t61.*t63.*t82.*4.0))-t103.*(l12y.*t3.*t4-l12z.*t4.*t12)+t9.*t19.*(l12z.*(t3.*t11-t2.*t12.*t13)+l12y.*t99)];
			mt2 = [t150.*(l12y.*(t52.*t53.*t54.^2.*t61.*t62.*4.0-t52.*t53.*t61.*t62.*t63.^2.*4.0)-l12x.*t11.*t13+l12z.*t3.*t4.*t11)+t103.*(l12x.*t4+l12z.*t3.*t13+l12y.*t12.*t13)+t9.*t19.*(-l12x.*t2.*t13+l12z.*t4.*t39+l12y.*t2.*t4.*t12),-t149.*(t27.*t118+t33.*t178)+t6.*t16.*(t33.*t118-t27.*t178),-t102.*t178+t27.*t149.*t177-t6.*t16.*t33.*t177,-t102.*(t166+t28.*(t141-t145))+t27.*t149.*t176-t6.*t16.*t33.*t176,-t102.*(t28.*t159+t34.*t161)-t27.*t149.*(t164-t167)+t6.*t16.*t33.*(t164-t167)];
			mt3 = [-t102.*(t28.*t151-t34.*t148)-t149.*(t33.*t119+t27.*t173)-t6.*t16.*(t27.*t119-t33.*t173),t149.*(t27.*(t156+t28.*(t133-t134))+t31.*t33.*t97)+t102.*(t28.*t154-t34.*(t133-t134))-t6.*t16.*(t33.*(t156+t28.*(t133-t134))-t27.*t31.*t97)];
			out1 = [mt1,mt2,mt3];
			end
			function out1 = ConsJac_3(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsJac_3
			%    OUT1 = ConsJac_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:22
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			q1dev__dt_0_ = in2(1,:);
			q1flex__dt_0_ = in2(2,:);
			q1sup__dt_0_ = in2(3,:);
			q2w1__dt_0_ = in2(4,:);
			q2w2__dt_0_ = in2(5,:);
			q2w3__dt_0_ = in2(6,:);
			q2w4__dt_0_ = in2(7,:);
			q2w5__dt_0_ = in2(8,:);
			q2w6__dt_0_ = in2(9,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf1Devx);
			t8 = cos(rf1Devy);
			t9 = sin(q1dev__dt_0_);
			t10 = sin(q1flex__dt_0_);
			t11 = sin(q1sup__dt_0_);
			t12 = sin(rf2Medx);
			t13 = sin(rf2Medy);
			t14 = sin(rf1Devx);
			t15 = sin(rf1Devy);
			t16 = l25+l13z;
			t17 = q2w1__dt_0_+riw1;
			t18 = q2w2__dt_0_+riw2;
			t19 = q2w3__dt_0_+riw3;
			t20 = q2w4__dt_0_+riw4;
			t21 = q2w5__dt_0_+riw5;
			t22 = q2w6__dt_0_+riw6;
			t37 = q1dev__dt_0_./2.0;
			t38 = q1flex__dt_0_./2.0;
			t39 = q1sup__dt_0_./2.0;
			t23 = cos(t17);
			t24 = cos(t18);
			t25 = cos(t19);
			t26 = cos(t20);
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = sin(t17);
			t30 = sin(t18);
			t31 = sin(t19);
			t32 = sin(t20);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = t2.*t3;
			t36 = t9.*t10;
			t40 = cos(t37);
			t41 = cos(t38);
			t42 = cos(t39);
			t43 = sin(t37);
			t44 = sin(t38);
			t45 = sin(t39);
			t46 = l13x.*t28;
			t47 = l13y.*t28;
			t48 = l13x.*t34;
			t49 = l13y.*t34;
			t50 = t11.*t36;
			t51 = t11.*t35;
			t52 = t16.*t27;
			t53 = t16.*t33;
			t55 = t40.^2;
			t56 = t41.^2;
			t54 = -t49;
			t57 = t47+t48;
			t58 = t35+t50;
			t59 = t36+t51;
			t60 = l23+t57;
			t61 = t46+t54;
			t70 = t26.*t33.*t57;
			t71 = t32.*t33.*t57;
			t62 = l24+t61;
			t63 = t26.*t60;
			t64 = t32.*t60;
			t65 = t26.*t61;
			t66 = t32.*t61;
			t73 = -t71;
			t67 = -t63;
			t68 = t27.*t62;
			t69 = t33.*t62;
			t84 = t66+t70;
			t85 = t65+t73;
			t72 = -t69;
			t74 = t53+t68;
			t86 = t25.*t84;
			t87 = t31.*t84;
			t88 = t25.*t85;
			t89 = t31.*t85;
			t75 = t52+t72;
			t78 = t25.*t26.*t74;
			t79 = t25.*t32.*t74;
			t80 = t26.*t31.*t74;
			t81 = t31.*t32.*t74;
			t90 = -t87;
			t106 = t86+t89;
			t109 = -t24.*(t87-t88);
			t76 = t26.*t75;
			t77 = t32.*t75;
			t83 = -t81;
			t102 = t79+t80;
			t107 = t88+t90;
			t108 = t30.*t106;
			t82 = -t77;
			t91 = t64+t76;
			t92 = t67+t77;
			t103 = t78+t83;
			t104 = t24.*t102;
			t110 = -t108;
			t93 = l22+t63+t82;
			t94 = t25.*t91;
			t95 = t31.*t91;
			t96 = -t25.*(t63+t82);
			t97 = -t31.*(t63+t82);
			t101 = t31.*(t63+t82);
			t105 = t30.*t103;
			t126 = t109+t110;
			t98 = -t95;
			t99 = t25.*t93;
			t100 = t31.*t93;
			t111 = t95+t96;
			t113 = t94+t101;
			t125 = t104+t105;
			t112 = t94+t100;
			t114 = t98+t99;
			t116 = t30.*t111;
			t119 = t24.*t113;
			t121 = -t30.*(t95-t99);
			t115 = l21+t114;
			t117 = t24.*t112;
			t118 = t30.*t112;
			t120 = -t116;
			t122 = t24.*t115;
			t123 = t30.*t115;
			t127 = t119+t120;
			t128 = t117+t121;
			t124 = -t122;
			t129 = t117+t123;
			t130 = t118+t124;
			mt1 = [t15.*(-l12z.*(t2.*t10-t3.*t9.*t11)+l12y.*t58+l12x.*t4.*t9)+t8.*t14.*(l12z.*t59+l12y.*(t40.*t43.*2.0-t40.*t43.*t56.*4.0-t41.*t42.*t43.^2.*t44.*t45.*4.0+t41.*t42.*t44.*t45.*t55.*4.0)+l12x.*t2.*t4),-t15.*(l12z.*(t3.*t9-t2.*t10.*t11)+l12y.*t59)+t7.*t8.*(l12y.*t3.*t4-l12z.*t4.*t10)-t8.*t14.*(l12z.*t58-l12y.*(t41.*t44.*2.0-t41.*t44.*t55.*4.0-t40.*t42.*t43.*t44.^2.*t45.*4.0+t40.*t42.*t43.*t45.*t56.*4.0))];
			mt2 = [-t15.*(-l12x.*t2.*t11+l12z.*t4.*t35+l12y.*t2.*t4.*t10)-t7.*t8.*(l12x.*t4+l12z.*t3.*t11+l12y.*t10.*t11)+t8.*t14.*(l12y.*(t40.*t41.*t42.^2.*t43.*t44.*4.0-t40.*t41.*t43.*t44.*t45.^2.*4.0)-l12x.*t9.*t11+l12z.*t3.*t4.*t9),-t13.*(t29.*t74-t23.*t130)-t6.*t12.*(t23.*t74+t29.*t130),t5.*t6.*t130+t13.*t29.*t129+t6.*t12.*t23.*t129,t5.*t6.*(t118+t24.*(t95-t99))+t13.*t29.*t128+t6.*t12.*t23.*t128,t5.*t6.*(t24.*t111+t30.*t113)-t13.*t29.*(t116-t119)-t6.*t12.*t23.*(t116-t119)];
			mt3 = [t13.*(t23.*t75-t29.*t125)-t6.*t12.*(t29.*t75+t23.*t125)+t5.*t6.*(t24.*t103-t30.*t102),t13.*(t29.*(t108+t24.*(t87-t88))-t23.*t27.*t57)-t5.*t6.*(t24.*t106-t30.*(t87-t88))+t6.*t12.*(t23.*(t108+t24.*(t87-t88))+t27.*t29.*t57)];
			out1 = [mt1,mt2,mt3];
			end
			function out1 = ConsJac_4(t,in2,in3,in4,in5,in6)
			%ConsJac_4
			%    OUT1 = ConsJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:26
			out1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			end
			function out1 = ConsJac_5(t,in2,in3,in4,in5,in6)
			%ConsJac_5
			%    OUT1 = ConsJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:27
			out1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			end
			function out1 = ConsJac_6(t,in2,in3,in4,in5,in6)
			%ConsJac_6
			%    OUT1 = ConsJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:28
			out1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
			end
			function out1 = ConsJac_7(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsJac_7
			%    OUT1 = ConsJac_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:29
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
			out1 = [t14.*t16.*t20.*(-1.0./2.0)+(t15.*t19.*t21)./2.0-(t17.*t24.*t54.*t57)./4.0+(t18.*t23.*t55.*t57)./4.0,t15.*t16.*t19.*(-1.0./2.0)+(t14.*t20.*t21)./2.0-(t18.*t23.*t54.*t57)./4.0+(t17.*t24.*t55.*t57)./4.0,t40.*(-1.0./2.0)+t45-(t17.*t24.*t57.*(t5./2.0+t6.*t22+t15.*t16.*t19.*2.0+t14.*t20.*t21.*(3.0./2.0)-t15.*t19.*t30+(t4.*t14.*t20.*t21)./2.0))./4.0-(t18.*t23.*t57.*(t6./2.0+t5.*t22+t14.*t16.*t20.*2.0+t15.*t19.*t21.*(3.0./2.0)-t14.*t20.*t30+(t4.*t15.*t19.*t21)./2.0))./4.0,0.0,0.0,0.0,0.0,0.0,0.0];
			end
			function out1 = ConsCor_1(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsCor_1
			%    OUT1 = ConsCor_1(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:15
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_0_ = in2(4,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_0_ = in2(5,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_0_ = in2(6,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_0_ = in2(7,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_0_ = in2(8,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_0_ = in2(9,:);
			q2w6__dt_1_ = in2(18,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = cos(rf1Devx);
			t9 = cos(rf1Devy);
			t10 = cos(rf1Devz);
			t11 = sin(q1dev__dt_0_);
			t12 = sin(q1flex__dt_0_);
			t13 = sin(q1sup__dt_0_);
			t14 = sin(rf2Medx);
			t15 = sin(rf2Medy);
			t16 = sin(rf2Medz);
			t17 = sin(rf1Devx);
			t18 = sin(rf1Devy);
			t19 = sin(rf1Devz);
			t20 = l25+l13z;
			t21 = q2w1__dt_0_+riw1;
			t22 = q2w2__dt_0_+riw2;
			t23 = q2w3__dt_0_+riw3;
			t24 = q2w4__dt_0_+riw4;
			t25 = q2w5__dt_0_+riw5;
			t26 = q2w6__dt_0_+riw6;
			t47 = q1dev__dt_0_./2.0;
			t48 = q1flex__dt_0_./2.0;
			t49 = q1sup__dt_0_./2.0;
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = cos(t23);
			t30 = cos(t24);
			t31 = cos(t25);
			t32 = cos(t26);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = sin(t23);
			t36 = sin(t24);
			t37 = sin(t25);
			t38 = sin(t26);
			t39 = t2.*t3;
			t40 = t2.*t12;
			t41 = t3.*t11;
			t42 = t5.*t16;
			t43 = t8.*t19;
			t44 = t11.*t12;
			t45 = t14.*t16;
			t46 = t17.*t19;
			t50 = cos(t47);
			t51 = cos(t48);
			t52 = cos(t49);
			t53 = sin(t47);
			t54 = sin(t48);
			t55 = sin(t49);
			t58 = l12x.*t2.*t4;
			t59 = l12z.*t3.*t4;
			t62 = l12x.*t2.*t13;
			t63 = l12x.*t4.*t11;
			t64 = l12y.*t3.*t13;
			t65 = l12y.*t4.*t12;
			t66 = l12x.*t11.*t13;
			t67 = l12z.*t12.*t13;
			t68 = t7.*t14.*t15;
			t69 = t10.*t17.*t18;
			t72 = t5.*t7.*t15;
			t73 = t8.*t10.*t18;
			t56 = l13x.*t32;
			t57 = l13y.*t32;
			t60 = l13x.*t38;
			t61 = l13y.*t38;
			t70 = t13.*t44;
			t71 = t13.*t39;
			t74 = t13.*t40;
			t75 = t13.*t41;
			t76 = t20.*t31;
			t77 = l12y.*t4.*t39;
			t78 = l12z.*t4.*t39;
			t79 = t20.*t37;
			t80 = l12z.*t4.*t40;
			t81 = l12z.*t4.*t41;
			t82 = l12y.*t4.*t44;
			t83 = l12z.*t4.*t44;
			t85 = t50.^2;
			t86 = t51.^2;
			t87 = t52.^2;
			t88 = -t66;
			t89 = -t67;
			t90 = t53.^2;
			t91 = t54.^2;
			t92 = t55.^2;
			t93 = -t68;
			t94 = -t69;
			t102 = t45+t72;
			t103 = t46+t73;
			t144 = t50.*t51.*t53.*t54.*4.0;
			t188 = t50.*t51.*t52.*t53.*t54.*t55.*8.0;
			t84 = -t61;
			t95 = -t74;
			t96 = -t75;
			t97 = -t78;
			t98 = -t80;
			t99 = t57+t60;
			t100 = t39+t70;
			t101 = t44+t71;
			t110 = t64+t89;
			t111 = t42+t93;
			t112 = t43+t94;
			t120 = t85.*t86.*2.0;
			t142 = t81+t82+t88;
			t151 = t50.*t53.*t86.*t92.*2.0;
			t152 = t50.*t53.*t87.*t91.*2.0;
			t153 = t51.*t54.*t87.*t90.*2.0;
			t154 = t51.*t54.*t85.*t92.*2.0;
			t155 = t52.*t55.*t85.*t91.*2.0;
			t156 = t52.*t55.*t86.*t90.*2.0;
			t157 = t50.*t53.*t91.*t92.*2.0;
			t158 = t51.*t54.*t90.*t92.*2.0;
			t159 = t52.*t55.*t90.*t91.*2.0;
			t162 = t50.*t53.*t86.*t87.*2.0;
			t163 = t51.*t54.*t85.*t87.*2.0;
			t104 = l23+t99;
			t105 = l12z.*t100;
			t106 = l12z.*t101;
			t107 = t56+t84;
			t108 = t40+t96;
			t109 = t41+t95;
			t118 = t30.*t99;
			t119 = t36.*t99;
			t121 = t77+t98;
			t129 = t27.*t37.*t99;
			t130 = t31.*t33.*t99;
			t132 = t33.*t37.*t99;
			t134 = t27.*t31.*t99;
			t150 = t103.*t110;
			t161 = t9.*t10.*t142;
			t164 = t52.*t55.*t120;
			t165 = -t155;
			t166 = -t156;
			t167 = -t157;
			t168 = -t158;
			t170 = -t162;
			t171 = -t163;
			t113 = l24+t107;
			t114 = l12y.*t108;
			t115 = l12y.*t109;
			t116 = l12z.*t108;
			t117 = l12z.*t109;
			t122 = t30.*t104;
			t123 = t36.*t104;
			t124 = t30.*t107;
			t125 = t36.*t107;
			t131 = t37.*t118;
			t133 = t37.*t119;
			t136 = t9.*t10.*t121;
			t137 = -t130;
			t139 = -t132;
			t145 = t29.*t31.*t118;
			t146 = t29.*t31.*t119;
			t147 = t31.*t35.*t118;
			t148 = t31.*t35.*t119;
			t169 = -t150;
			t173 = -t161;
			t253 = t151+t152+t167+t170;
			t254 = t153+t154+t168+t171;
			t282 = t144+t159+t164+t165+t166;
			t126 = -t122;
			t127 = t31.*t113;
			t128 = t37.*t113;
			t138 = t37.*t124;
			t140 = t37.*t125;
			t141 = -t133;
			t149 = -t148;
			t174 = t105+t114;
			t197 = t125+t131;
			t210 = t146+t147;
			t257 = l12y.*t253;
			t258 = l12y.*t254;
			t284 = l12y.*t282;
			t135 = -t128;
			t143 = -t138;
			t160 = t79+t127;
			t181 = t9.*t10.*t174;
			t196 = t118+t140;
			t198 = t124+t141;
			t201 = t29.*t197;
			t203 = t35.*t197;
			t211 = t145+t149;
			t213 = t28.*t210;
			t214 = t34.*t210;
			t268 = t83+t257;
			t279 = t62+t97+t258;
			t286 = t112.*(t62-t78+t258);
			t288 = t117+t284;
			t172 = t76+t135;
			t175 = t27.*t160;
			t176 = t33.*t160;
			t182 = t29.*t30.*t160;
			t183 = t29.*t36.*t160;
			t184 = t30.*t35.*t160;
			t185 = t35.*t36.*t160;
			t189 = -t181;
			t199 = t119+t143;
			t200 = t29.*t196;
			t202 = t35.*t196;
			t204 = t29.*t198;
			t206 = t35.*t198;
			t208 = -t203;
			t217 = t28.*t211;
			t218 = t34.*t211;
			t219 = -t214;
			t283 = t112.*t268;
			t295 = t112.*t288;
			t298 = t173+t286;
			t177 = t27.*t172;
			t178 = t30.*t172;
			t179 = t33.*t172;
			t180 = t36.*t172;
			t194 = -t185;
			t205 = t29.*t199;
			t207 = t35.*t199;
			t228 = t183+t184;
			t239 = t201+t206;
			t242 = t204+t208;
			t247 = -t28.*(t203-t204);
			t248 = -t34.*(t203-t204);
			t251 = t213+t218;
			t252 = t217+t219;
			t275 = -t102.*(t214-t217);
			t300 = t136+t169+t283;
			t311 = t189+t295;
			t186 = -t177;
			t187 = -t180;
			t190 = t29.*t178;
			t191 = t29.*t180;
			t192 = t35.*t178;
			t193 = t35.*t180;
			t209 = -t207;
			t212 = t123+t178;
			t215 = t126+t180;
			t229 = t182+t194;
			t231 = t28.*t228;
			t232 = t34.*t228;
			t240 = t202+t205;
			t243 = t28.*t239;
			t244 = t34.*t239;
			t255 = t27.*t251;
			t256 = t33.*t251;
			t195 = t35.*t187;
			t216 = l22+t122+t187;
			t220 = t29.*t212;
			t221 = t35.*t212;
			t222 = -t29.*(t122+t187);
			t223 = -t35.*(t122+t187);
			t227 = t35.*(t122+t187);
			t230 = t191+t192;
			t234 = t28.*t229;
			t235 = t34.*t229;
			t236 = -t232;
			t241 = t200+t209;
			t245 = t34.*t240;
			t249 = -t244;
			t280 = t129+t256;
			t281 = t139+t255;
			t287 = -t111.*(t132-t255);
			t301 = t243+t248;
			t306 = -t27.*(t244+t28.*(t203-t204));
			t307 = -t33.*(t244+t28.*(t203-t204));
			t312 = -t102.*(t244+t28.*(t203-t204));
			t317 = -t6.*t7.*(t130+t27.*(t244+t28.*(t203-t204)));
			t318 = t6.*t7.*(t130+t27.*(t244+t28.*(t203-t204)));
			t224 = -t221;
			t225 = t29.*t216;
			t226 = t35.*t216;
			t233 = t190+t195;
			t237 = t28.*t230;
			t246 = t28.*t241;
			t250 = -t245;
			t259 = t221+t222;
			t261 = t220+t227;
			t285 = t6.*t7.*t280;
			t289 = t231+t235;
			t290 = t234+t236;
			t294 = -t6.*t7.*t33.*(t232-t234);
			t299 = -t27.*t111.*(t232-t234);
			t303 = t247+t249;
			t309 = t6.*t7.*t33.*t301;
			t313 = t27.*t111.*t301;
			t315 = t134+t307;
			t316 = t137+t306;
			t238 = t34.*t233;
			t260 = t220+t226;
			t262 = t224+t225;
			t264 = t28.*t259;
			t265 = t34.*t259;
			t269 = t28.*t261;
			t270 = t34.*t261;
			t272 = -t28.*(t221-t225);
			t273 = -t34.*(t221-t225);
			t291 = t33.*t289;
			t292 = t27.*t289;
			t296 = t102.*t289;
			t302 = t246+t250;
			t319 = t111.*t315;
			t342 = t275+t285+t287;
			t350 = t309+t312+t313;
			t263 = l21+t262;
			t266 = t28.*t260;
			t267 = t34.*t260;
			t271 = -t265;
			t293 = t237+t238;
			t297 = -t296;
			t304 = t179+t292;
			t305 = t186+t291;
			t314 = -t111.*(t177-t291);
			t320 = t264+t270;
			t328 = -t6.*t7.*t27.*(t265-t269);
			t332 = t6.*t7.*t27.*(t265-t269);
			t335 = -t102.*(t265-t269);
			t339 = -t33.*t111.*(t265-t269);
			t344 = t318+t319;
			t346 = q2w5__dt_1_.*(t296+t27.*t111.*(t232-t234)+t6.*t7.*t33.*(t232-t234)).*-2.0;
			t351 = q2w6__dt_1_.*t350.*2.0;
			t274 = -t267;
			t276 = t28.*t263;
			t277 = t34.*t263;
			t308 = t6.*t7.*t304;
			t321 = t269+t271;
			t322 = t266+t273;
			t326 = t6.*t7.*t33.*t320;
			t331 = -t6.*t7.*t33.*(t267+t28.*(t221-t225));
			t337 = t27.*t111.*t320;
			t341 = -t27.*t111.*(t267+t28.*(t221-t225));
			t345 = t294+t297+t299;
			t348 = t332+t339;
			t352 = -t351;
			t278 = -t276;
			t310 = -t308;
			t323 = t266+t277;
			t324 = t272+t274;
			t327 = t6.*t7.*t27.*t322;
			t334 = t102.*t322;
			t338 = t33.*t111.*t322;
			t353 = t326+t335+t337;
			t325 = t267+t278;
			t329 = t6.*t7.*t27.*t323;
			t330 = -t327;
			t336 = -t334;
			t340 = t33.*t111.*t323;
			t343 = t310+t314;
			t354 = q2w4__dt_1_.*t353.*2.0;
			t356 = q2w3__dt_1_.*(t334+t27.*t111.*(t267+t28.*(t221-t225))+t6.*t7.*t33.*(t267+t28.*(t221-t225))).*-2.0;
			t357 = q2w3__dt_1_.*(t334+t27.*t111.*(t267+t28.*(t221-t225))+t6.*t7.*t33.*(t267+t28.*(t221-t225))).*2.0;
			t333 = -t329;
			t347 = t330+t338;
			t355 = t331+t336+t341;
			et1 = -q1sup__dt_1_.*(q1flex__dt_1_.*t300.*-2.0+q1sup__dt_1_.*(-t112.*(t63+l12z.*t75+l12y.*t188)+t103.*(t59+t65-l12x.*t13)+t9.*t10.*(t58+l12z.*t71+l12y.*t74)).*2.0+q1dev__dt_1_.*(t161-t286).*2.0);
			et2 = q2w5__dt_1_.*(q2w1__dt_1_.*(t308+t111.*(t177-t291)).*2.0+q2w5__dt_1_.*(t102.*(t28.*t233-t34.*t230)-t111.*(t176-t27.*t293)+t6.*t7.*(t175+t33.*t293)).*2.0-q2w2__dt_1_.*(t296+t27.*t111.*(t232-t234)+t6.*t7.*t33.*(t232-t234)).*2.0-q2w3__dt_1_.*(t296+t27.*t111.*(t232-t234)+t6.*t7.*t33.*(t232-t234)).*2.0-q2w4__dt_1_.*(t296+t27.*t111.*(t232-t234)+t6.*t7.*t33.*(t232-t234)).*2.0+q2w6__dt_1_.*(-t285+t111.*(t132-t255)+t102.*(t214-t217)).*2.0);
			et3 = q1flex__dt_1_.*(q1sup__dt_1_.*t300.*2.0+q1dev__dt_1_.*(t181-t295).*2.0-q1flex__dt_1_.*(t103.*(t59+t65)+t112.*(t116-l12y.*(-t86+t91+t120+t188-t85.*t91.*2.0))+t9.*t10.*(t106-t115)).*2.0)-q2w6__dt_1_.*(q2w6__dt_1_.*(t111.*(t27.*(t245-t246)-t31.*t33.*t107)+t102.*(t28.*t240+t34.*t241)+t6.*t7.*(t33.*(t245-t246)+t27.*t31.*t107)).*-2.0+q2w1__dt_1_.*t344.*2.0+q2w2__dt_1_.*t350.*2.0+q2w3__dt_1_.*t350.*2.0+q2w4__dt_1_.*t350.*2.0-q2w5__dt_1_.*(-t285+t111.*(t132-t255)+t102.*(t214-t217)).*2.0);
			et4 = -q2w4__dt_1_.*(t351-t354-q2w2__dt_1_.*t353.*2.0-q2w3__dt_1_.*t353.*2.0+q2w1__dt_1_.*(t328+t33.*t111.*(t265-t269)).*2.0+q2w5__dt_1_.*(t296+t27.*t111.*(t232-t234)+t6.*t7.*t33.*(t232-t234)).*2.0)-q2w1__dt_1_.*(q2w5__dt_1_.*(t308+t111.*(t177-t291)).*-2.0+q2w6__dt_1_.*t344.*2.0+q2w2__dt_1_.*(t329-t340).*2.0+q2w3__dt_1_.*(t327-t338).*2.0+q2w1__dt_1_.*(t111.*(t176-t27.*t325)-t6.*t7.*(t175+t33.*t325)).*2.0+q2w4__dt_1_.*(t328+t33.*t111.*(t265-t269)).*2.0)+q2w3__dt_1_.*(t346+t352+t354+t357+q2w2__dt_1_.*(t334+t27.*t111.*(t267+t28.*(t221-t225))+t6.*t7.*t33.*(t267+t28.*(t221-t225))).*2.0-q2w1__dt_1_.*(t327-t338).*2.0);
			et5 = q1dev__dt_1_.*(q1flex__dt_1_.*(t181-t295).*2.0-q1sup__dt_1_.*(t161-t286).*2.0+q1dev__dt_1_.*(t112.*(t63-t116+l12y.*(-t85+t90+t120+t188-t86.*t90.*2.0))-t9.*t10.*(t58+t106-t115)).*2.0)+q2w2__dt_1_.*(t346+t352+t354+t357+q2w2__dt_1_.*(t102.*t323+t27.*t111.*t325+t6.*t7.*t33.*t325).*2.0-q2w1__dt_1_.*(t329-t340).*2.0);
			out1 = et1+et2+et3+et4+et5;
			end
			function out1 = ConsCor_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsCor_2
			%    OUT1 = ConsCor_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:19
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_0_ = in2(4,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_0_ = in2(5,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_0_ = in2(6,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_0_ = in2(7,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_0_ = in2(8,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_0_ = in2(9,:);
			q2w6__dt_1_ = in2(18,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = cos(rf1Devx);
			t9 = cos(rf1Devy);
			t10 = cos(rf1Devz);
			t11 = sin(q1dev__dt_0_);
			t12 = sin(q1flex__dt_0_);
			t13 = sin(q1sup__dt_0_);
			t14 = sin(rf2Medx);
			t15 = sin(rf2Medy);
			t16 = sin(rf2Medz);
			t17 = sin(rf1Devx);
			t18 = sin(rf1Devy);
			t19 = sin(rf1Devz);
			t20 = l25+l13z;
			t21 = q2w1__dt_0_+riw1;
			t22 = q2w2__dt_0_+riw2;
			t23 = q2w3__dt_0_+riw3;
			t24 = q2w4__dt_0_+riw4;
			t25 = q2w5__dt_0_+riw5;
			t26 = q2w6__dt_0_+riw6;
			t45 = q1dev__dt_0_./2.0;
			t46 = q1flex__dt_0_./2.0;
			t47 = q1sup__dt_0_./2.0;
			t48 = rf2Medx./2.0;
			t49 = rf2Medy./2.0;
			t50 = rf2Medz./2.0;
			t51 = rf1Devx./2.0;
			t52 = rf1Devy./2.0;
			t53 = rf1Devz./2.0;
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = cos(t23);
			t30 = cos(t24);
			t31 = cos(t25);
			t32 = cos(t26);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = sin(t23);
			t36 = sin(t24);
			t37 = sin(t25);
			t38 = sin(t26);
			t39 = t2.*t3;
			t40 = t2.*t12;
			t41 = t3.*t11;
			t42 = t7.*t14;
			t43 = t10.*t17;
			t44 = t11.*t12;
			t54 = cos(t45);
			t55 = cos(t46);
			t56 = cos(t47);
			t57 = cos(t48);
			t58 = cos(t49);
			t59 = cos(t50);
			t60 = cos(t51);
			t61 = cos(t52);
			t62 = cos(t53);
			t63 = sin(t45);
			t64 = sin(t46);
			t65 = sin(t47);
			t66 = sin(t48);
			t67 = sin(t49);
			t68 = sin(t50);
			t69 = sin(t51);
			t70 = sin(t52);
			t71 = sin(t53);
			t74 = l12x.*t2.*t4;
			t75 = l12z.*t3.*t4;
			t78 = l12x.*t2.*t13;
			t79 = l12x.*t4.*t11;
			t80 = l12y.*t3.*t13;
			t81 = l12y.*t4.*t12;
			t82 = l12x.*t11.*t13;
			t83 = l12z.*t12.*t13;
			t84 = t5.*t15.*t16;
			t85 = t8.*t18.*t19;
			t72 = l13x.*t32;
			t73 = l13y.*t32;
			t76 = l13x.*t38;
			t77 = l13y.*t38;
			t86 = t13.*t44;
			t87 = t13.*t39;
			t88 = t13.*t40;
			t89 = t13.*t41;
			t90 = t20.*t31;
			t91 = l12y.*t4.*t39;
			t92 = l12z.*t4.*t39;
			t93 = t20.*t37;
			t94 = l12z.*t4.*t40;
			t95 = l12z.*t4.*t41;
			t96 = l12y.*t4.*t44;
			t97 = l12z.*t4.*t44;
			t99 = t54.^2;
			t100 = t55.^2;
			t101 = t56.^2;
			t102 = -t82;
			t103 = -t83;
			t104 = t57.^2;
			t105 = t59.^2;
			t106 = t60.^2;
			t107 = t62.^2;
			t108 = t63.^2;
			t109 = t64.^2;
			t110 = t65.^2;
			t111 = -t84;
			t112 = -t85;
			t170 = t54.*t55.*t63.*t64.*4.0;
			t212 = t54.*t55.*t56.*t63.*t64.*t65.*8.0;
			t213 = t57.*t58.*t59.*t66.*t67.*t68.*8.0;
			t214 = t60.*t61.*t62.*t69.*t70.*t71.*8.0;
			t98 = -t77;
			t113 = -t88;
			t114 = -t89;
			t115 = -t92;
			t116 = -t94;
			t117 = t104.*2.0;
			t118 = t105.*2.0;
			t119 = t106.*2.0;
			t120 = t107.*2.0;
			t125 = t73+t76;
			t126 = t39+t86;
			t127 = t44+t87;
			t134 = t80+t103;
			t135 = t42+t111;
			t136 = t43+t112;
			t144 = t99.*t100.*2.0;
			t145 = t104.*t105.*4.0;
			t146 = t106.*t107.*4.0;
			t168 = t95+t96+t102;
			t176 = t54.*t63.*t100.*t110.*2.0;
			t177 = t54.*t63.*t101.*t109.*2.0;
			t178 = t55.*t64.*t101.*t108.*2.0;
			t179 = t55.*t64.*t99.*t110.*2.0;
			t180 = t56.*t65.*t99.*t109.*2.0;
			t181 = t56.*t65.*t100.*t108.*2.0;
			t182 = t54.*t63.*t109.*t110.*2.0;
			t183 = t55.*t64.*t108.*t110.*2.0;
			t184 = t56.*t65.*t108.*t109.*2.0;
			t187 = t54.*t63.*t100.*t101.*2.0;
			t188 = t55.*t64.*t99.*t101.*2.0;
			t121 = -t117;
			t122 = -t118;
			t123 = -t119;
			t124 = -t120;
			t128 = l23+t125;
			t129 = l12z.*t126;
			t130 = l12z.*t127;
			t131 = t72+t98;
			t132 = t40+t114;
			t133 = t41+t113;
			t142 = t30.*t125;
			t143 = t36.*t125;
			t147 = t91+t116;
			t155 = t27.*t37.*t125;
			t156 = t31.*t33.*t125;
			t158 = t33.*t37.*t125;
			t160 = t27.*t31.*t125;
			t186 = t9.*t19.*t168;
			t189 = t56.*t65.*t144;
			t190 = -t180;
			t191 = -t181;
			t192 = -t182;
			t193 = -t183;
			t194 = t134.*t136;
			t195 = -t187;
			t196 = -t188;
			t137 = l24+t131;
			t138 = l12y.*t132;
			t139 = l12y.*t133;
			t140 = l12z.*t132;
			t141 = l12z.*t133;
			t148 = t30.*t128;
			t149 = t36.*t128;
			t150 = t30.*t131;
			t151 = t36.*t131;
			t157 = t37.*t142;
			t159 = t37.*t143;
			t162 = -t156;
			t164 = t9.*t19.*t147;
			t165 = -t158;
			t171 = t29.*t31.*t142;
			t172 = t29.*t31.*t143;
			t173 = t31.*t35.*t142;
			t174 = t31.*t35.*t143;
			t254 = t121+t122+t145+t213+1.0;
			t255 = t123+t124+t146+t214+1.0;
			t280 = t176+t177+t192+t195;
			t281 = t178+t179+t193+t196;
			t309 = t170+t184+t189+t190+t191;
			t152 = -t148;
			t153 = t31.*t137;
			t154 = t37.*t137;
			t163 = t37.*t150;
			t166 = t37.*t151;
			t167 = -t159;
			t175 = -t174;
			t198 = t129+t138;
			t222 = t151+t157;
			t235 = t172+t173;
			t284 = l12y.*t280;
			t285 = l12y.*t281;
			t310 = l12y.*t309;
			t161 = -t154;
			t169 = -t163;
			t185 = t93+t153;
			t205 = t9.*t19.*t198;
			t221 = t142+t166;
			t223 = t150+t167;
			t226 = t29.*t222;
			t228 = t35.*t222;
			t236 = t171+t175;
			t238 = t28.*t235;
			t239 = t34.*t235;
			t295 = t97+t284;
			t306 = t78+t115+t285;
			t313 = t141+t310;
			t333 = t255.*(t78-t92+t285);
			t197 = t90+t161;
			t199 = t27.*t185;
			t200 = t33.*t185;
			t206 = t29.*t30.*t185;
			t207 = t29.*t36.*t185;
			t208 = t30.*t35.*t185;
			t209 = t35.*t36.*t185;
			t224 = t143+t169;
			t225 = t29.*t221;
			t227 = t35.*t221;
			t229 = t29.*t223;
			t231 = t35.*t223;
			t233 = -t228;
			t242 = t28.*t236;
			t243 = t34.*t236;
			t244 = -t239;
			t325 = t255.*t295;
			t339 = t255.*t313;
			t340 = t186+t333;
			t201 = t27.*t197;
			t202 = t30.*t197;
			t203 = t33.*t197;
			t204 = t36.*t197;
			t219 = -t209;
			t230 = t29.*t224;
			t232 = t35.*t224;
			t253 = t207+t208;
			t266 = t226+t231;
			t269 = t229+t233;
			t274 = -t28.*(t228-t229);
			t275 = -t34.*(t228-t229);
			t278 = t238+t243;
			t279 = t242+t244;
			t304 = -t135.*(t239-t242);
			t327 = -t325;
			t358 = t205+t339;
			t210 = -t201;
			t211 = -t204;
			t215 = t29.*t202;
			t216 = t29.*t204;
			t217 = t35.*t202;
			t218 = t35.*t204;
			t234 = -t232;
			t237 = t149+t202;
			t240 = t152+t204;
			t256 = t206+t219;
			t258 = t28.*t253;
			t259 = t34.*t253;
			t267 = t227+t230;
			t270 = t28.*t266;
			t271 = t34.*t266;
			t282 = t27.*t278;
			t283 = t33.*t278;
			t349 = t164+t194+t327;
			t220 = t35.*t211;
			t241 = l22+t148+t211;
			t245 = t29.*t237;
			t246 = t35.*t237;
			t247 = -t29.*(t148+t211);
			t248 = -t35.*(t148+t211);
			t252 = t35.*(t148+t211);
			t257 = t216+t217;
			t261 = t28.*t256;
			t262 = t34.*t256;
			t263 = -t259;
			t268 = t225+t234;
			t272 = t34.*t267;
			t276 = -t271;
			t307 = t155+t283;
			t308 = t165+t282;
			t321 = t270+t275;
			t328 = -t27.*(t271+t28.*(t228-t229));
			t329 = -t33.*(t271+t28.*(t228-t229));
			t334 = -t135.*(t271+t28.*(t228-t229));
			t335 = -t254.*(t158-t282);
			t338 = -t6.*t16.*(t156+t27.*(t271+t28.*(t228-t229)));
			t249 = -t246;
			t250 = t29.*t241;
			t251 = t35.*t241;
			t260 = t215+t220;
			t264 = t28.*t257;
			t273 = t28.*t268;
			t277 = -t272;
			t286 = t246+t247;
			t288 = t245+t252;
			t311 = t6.*t16.*t307;
			t314 = t258+t262;
			t315 = t261+t263;
			t319 = -t6.*t16.*t33.*(t259-t261);
			t323 = t274+t276;
			t331 = t6.*t16.*t33.*t321;
			t336 = t160+t329;
			t337 = t162+t328;
			t341 = -t27.*t254.*(t259-t261);
			t342 = t27.*t254.*(t259-t261);
			t359 = t27.*t254.*t321;
			t265 = t34.*t260;
			t287 = t245+t251;
			t289 = t249+t250;
			t291 = t28.*t286;
			t292 = t34.*t286;
			t296 = t28.*t288;
			t297 = t34.*t288;
			t299 = -t28.*(t246-t250);
			t300 = -t34.*(t246-t250);
			t312 = -t311;
			t316 = t33.*t314;
			t317 = t27.*t314;
			t320 = t135.*t314;
			t322 = t273+t277;
			t332 = -t331;
			t361 = t254.*t336;
			t378 = q2w6__dt_1_.*(t331-t359+t135.*(t271+t28.*(t228-t229))).*-2.0;
			t379 = q2w6__dt_1_.*(t331-t359+t135.*(t271+t28.*(t228-t229))).*2.0;
			t290 = l21+t289;
			t293 = t28.*t287;
			t294 = t34.*t287;
			t298 = -t292;
			t318 = t264+t265;
			t324 = t203+t317;
			t326 = t210+t316;
			t343 = t291+t297;
			t352 = -t6.*t16.*t27.*(t292-t296);
			t357 = -t135.*(t292-t296);
			t360 = -t254.*(t201-t316);
			t364 = -t33.*t254.*(t292-t296);
			t369 = t304+t312+t335;
			t370 = t338+t361;
			t371 = t319+t320+t342;
			t377 = t332+t334+t359;
			t301 = -t294;
			t302 = t28.*t290;
			t303 = t34.*t290;
			t330 = t6.*t16.*t324;
			t344 = t296+t298;
			t345 = t293+t300;
			t350 = t6.*t16.*t33.*t343;
			t355 = -t6.*t16.*t33.*(t294+t28.*(t246-t250));
			t362 = t27.*t254.*t343;
			t366 = -t27.*t254.*(t294+t28.*(t246-t250));
			t367 = t27.*t254.*(t294+t28.*(t246-t250));
			t372 = t352+t364;
			t373 = q2w5__dt_1_.*t371.*2.0;
			t305 = -t302;
			t346 = t293+t303;
			t347 = t299+t301;
			t351 = t6.*t16.*t27.*t345;
			t353 = -t350;
			t356 = t135.*t345;
			t363 = t33.*t254.*t345;
			t368 = t330+t360;
			t375 = -t373;
			t381 = q2w4__dt_1_.*(t350-t362+t135.*(t292-t296)).*-2.0;
			t348 = t294+t305;
			t354 = t6.*t16.*t27.*t346;
			t365 = t33.*t254.*t346;
			t374 = t351+t363;
			t380 = t353+t357+t362;
			t382 = t355+t356+t367;
			t376 = t354+t365;
			t383 = q2w3__dt_1_.*t382.*2.0;
			et1 = -q1sup__dt_1_.*(q1dev__dt_1_.*t340.*2.0-q1flex__dt_1_.*t349.*2.0+q1sup__dt_1_.*(t255.*(t79+l12z.*t89+l12y.*t212)-t136.*(t75+t81-l12x.*t13)+t9.*t19.*(t74+l12y.*t88+l12z.*t87)).*2.0)-q2w6__dt_1_.*(q2w6__dt_1_.*(t254.*(t27.*(t272-t273)-t31.*t33.*t131)+t135.*(t28.*t267+t34.*t268)-t6.*t16.*(t33.*(t272-t273)+t27.*t31.*t131)).*2.0-q2w1__dt_1_.*t370.*2.0+q2w2__dt_1_.*(t331-t359+t135.*(t271+t28.*(t228-t229))).*2.0+q2w3__dt_1_.*(t331-t359+t135.*(t271+t28.*(t228-t229))).*2.0+q2w4__dt_1_.*(t331-t359+t135.*(t271+t28.*(t228-t229))).*2.0+q2w5__dt_1_.*(t311+t135.*(t239-t242)+t254.*(t158-t282)).*2.0);
			et2 = q2w5__dt_1_.*(q2w5__dt_1_.*(-t135.*(t28.*t260-t34.*t257)+t254.*(t200-t27.*t318)+t6.*t16.*(t199+t33.*t318)).*2.0+q2w1__dt_1_.*t368.*2.0+q2w2__dt_1_.*t371.*2.0+q2w3__dt_1_.*t371.*2.0+q2w4__dt_1_.*t371.*2.0-q2w6__dt_1_.*(t311+t135.*(t239-t242)+t254.*(t158-t282)).*2.0)+q1flex__dt_1_.*(q1dev__dt_1_.*t358.*2.0+q1sup__dt_1_.*t349.*2.0+q1flex__dt_1_.*(t136.*(t75+t81)+t255.*(t140-l12y.*(-t100+t109+t144+t212-t99.*t109.*2.0))-t9.*t19.*(t130-t139)).*2.0)+q2w1__dt_1_.*(q2w2__dt_1_.*t376.*-2.0-q2w3__dt_1_.*t374.*2.0+q2w5__dt_1_.*t368.*2.0+q2w6__dt_1_.*t370.*2.0+q2w1__dt_1_.*(t254.*(t200-t27.*t348)+t6.*t16.*(t199+t33.*t348)).*2.0+q2w4__dt_1_.*(t33.*t254.*(t292-t296)+t6.*t16.*t27.*(t292-t296)).*2.0);
			et3 = -q1dev__dt_1_.*(q1flex__dt_1_.*t358.*-2.0+q1sup__dt_1_.*t340.*2.0+q1dev__dt_1_.*(t255.*(t79-t140+l12y.*(-t99+t108+t144+t212-t100.*t108.*2.0))+t9.*t19.*(t74+t130-t139)).*2.0)+q2w4__dt_1_.*(t373+t378+q2w1__dt_1_.*(t33.*t254.*(t292-t296)+t6.*t16.*t27.*(t292-t296)).*2.0+q2w2__dt_1_.*(t350-t362+t135.*(t292-t296)).*2.0+q2w3__dt_1_.*(t350-t362+t135.*(t292-t296)).*2.0+q2w4__dt_1_.*(t350-t362+t135.*(t292-t296)).*2.0)-q2w3__dt_1_.*(t375+t379+t381+t383+q2w1__dt_1_.*t374.*2.0+q2w2__dt_1_.*t382.*2.0)-q2w2__dt_1_.*(t375+t379+t381+t383+q2w1__dt_1_.*t376.*2.0+q2w2__dt_1_.*(t135.*t346+t27.*t254.*t348-t6.*t16.*t33.*t348).*2.0);
			out1 = et1+et2+et3;
			end
			function out1 = ConsCor_3(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsCor_3
			%    OUT1 = ConsCor_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:24
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_0_ = in2(4,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_0_ = in2(5,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_0_ = in2(6,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_0_ = in2(7,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_0_ = in2(8,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_0_ = in2(9,:);
			q2w6__dt_1_ = in2(18,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf1Devx);
			t8 = cos(rf1Devy);
			t9 = sin(q1dev__dt_0_);
			t10 = sin(q1flex__dt_0_);
			t11 = sin(q1sup__dt_0_);
			t12 = sin(rf2Medx);
			t13 = sin(rf2Medy);
			t14 = sin(rf1Devx);
			t15 = sin(rf1Devy);
			t16 = l25+l13z;
			t17 = q2w1__dt_0_+riw1;
			t18 = q2w2__dt_0_+riw2;
			t19 = q2w3__dt_0_+riw3;
			t20 = q2w4__dt_0_+riw4;
			t21 = q2w5__dt_0_+riw5;
			t22 = q2w6__dt_0_+riw6;
			t39 = q1dev__dt_0_./2.0;
			t40 = q1flex__dt_0_./2.0;
			t41 = q1sup__dt_0_./2.0;
			t23 = cos(t17);
			t24 = cos(t18);
			t25 = cos(t19);
			t26 = cos(t20);
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = sin(t17);
			t30 = sin(t18);
			t31 = sin(t19);
			t32 = sin(t20);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = t2.*t3;
			t36 = t2.*t10;
			t37 = t3.*t9;
			t38 = t9.*t10;
			t42 = cos(t39);
			t43 = cos(t40);
			t44 = cos(t41);
			t45 = sin(t39);
			t46 = sin(t40);
			t47 = sin(t41);
			t50 = l12x.*t2.*t4;
			t51 = l12z.*t3.*t4;
			t54 = l12x.*t2.*t11;
			t55 = l12x.*t4.*t9;
			t56 = l12y.*t3.*t11;
			t57 = l12y.*t4.*t10;
			t58 = l12x.*t9.*t11;
			t59 = l12z.*t10.*t11;
			t48 = l13x.*t28;
			t49 = l13y.*t28;
			t52 = l13x.*t34;
			t53 = l13y.*t34;
			t60 = t11.*t38;
			t61 = t11.*t35;
			t62 = t11.*t36;
			t63 = t11.*t37;
			t64 = t16.*t27;
			t65 = l12y.*t4.*t35;
			t66 = l12z.*t4.*t35;
			t67 = t16.*t33;
			t68 = l12z.*t4.*t36;
			t69 = l12z.*t4.*t37;
			t70 = l12y.*t4.*t38;
			t71 = l12z.*t4.*t38;
			t73 = t42.^2;
			t74 = t43.^2;
			t75 = t44.^2;
			t76 = -t58;
			t77 = -t59;
			t78 = t45.^2;
			t79 = t46.^2;
			t80 = t47.^2;
			t127 = t42.*t43.*t45.*t46.*4.0;
			t168 = t42.*t43.*t44.*t45.*t46.*t47.*8.0;
			t72 = -t53;
			t81 = -t62;
			t82 = -t63;
			t83 = -t66;
			t84 = -t68;
			t85 = t49+t52;
			t86 = t35+t60;
			t87 = t38+t61;
			t94 = t56+t77;
			t102 = t73.*t74.*2.0;
			t125 = t69+t70+t76;
			t134 = t42.*t45.*t74.*t80.*2.0;
			t135 = t42.*t45.*t75.*t79.*2.0;
			t136 = t43.*t46.*t75.*t78.*2.0;
			t137 = t43.*t46.*t73.*t80.*2.0;
			t138 = t44.*t47.*t73.*t79.*2.0;
			t139 = t44.*t47.*t74.*t78.*2.0;
			t140 = t42.*t45.*t79.*t80.*2.0;
			t141 = t43.*t46.*t78.*t80.*2.0;
			t142 = t44.*t47.*t78.*t79.*2.0;
			t144 = t42.*t45.*t74.*t75.*2.0;
			t145 = t43.*t46.*t73.*t75.*2.0;
			t88 = l23+t85;
			t89 = l12z.*t86;
			t90 = l12z.*t87;
			t91 = t48+t72;
			t92 = t36+t82;
			t93 = t37+t81;
			t100 = t26.*t85;
			t101 = t32.*t85;
			t103 = t65+t84;
			t108 = t7.*t8.*t94;
			t112 = t23.*t33.*t85;
			t113 = t27.*t29.*t85;
			t115 = t29.*t33.*t85;
			t118 = t23.*t27.*t85;
			t133 = t15.*t125;
			t146 = t44.*t47.*t102;
			t147 = -t138;
			t148 = -t139;
			t149 = -t140;
			t150 = -t141;
			t151 = -t144;
			t152 = -t145;
			t95 = l24+t91;
			t96 = l12y.*t92;
			t97 = l12y.*t93;
			t98 = l12z.*t92;
			t99 = l12z.*t93;
			t104 = t26.*t88;
			t105 = t32.*t88;
			t106 = t26.*t91;
			t107 = t32.*t91;
			t114 = t33.*t100;
			t116 = t33.*t101;
			t117 = t15.*t103;
			t120 = -t113;
			t122 = -t115;
			t128 = t25.*t27.*t100;
			t129 = t25.*t27.*t101;
			t130 = t27.*t31.*t100;
			t131 = t27.*t31.*t101;
			t232 = t134+t135+t149+t151;
			t233 = t136+t137+t150+t152;
			t265 = t127+t142+t146+t147+t148;
			t109 = -t104;
			t110 = t27.*t95;
			t111 = t33.*t95;
			t121 = t33.*t106;
			t123 = t33.*t107;
			t124 = -t116;
			t132 = -t131;
			t154 = t89+t96;
			t176 = t107+t114;
			t189 = t129+t130;
			t236 = l12y.*t232;
			t237 = l12y.*t233;
			t267 = l12y.*t265;
			t119 = -t111;
			t126 = -t121;
			t143 = t67+t110;
			t157 = t15.*t154;
			t175 = t100+t123;
			t177 = t106+t124;
			t180 = t25.*t176;
			t182 = t31.*t176;
			t190 = t128+t132;
			t192 = t24.*t189;
			t193 = t30.*t189;
			t249 = t71+t236;
			t260 = t54+t83+t237;
			t262 = t8.*t14.*(t54-t66+t237);
			t269 = t99+t267;
			t153 = t64+t119;
			t155 = t23.*t143;
			t156 = t29.*t143;
			t162 = t25.*t26.*t143;
			t163 = t25.*t32.*t143;
			t164 = t26.*t31.*t143;
			t165 = t31.*t32.*t143;
			t178 = t101+t126;
			t179 = t25.*t175;
			t181 = t31.*t175;
			t183 = t25.*t177;
			t185 = t31.*t177;
			t187 = -t182;
			t196 = t24.*t190;
			t197 = t30.*t190;
			t198 = -t193;
			t259 = t8.*t14.*t249;
			t264 = -t262;
			t272 = t8.*t14.*t269;
			t158 = t23.*t153;
			t159 = t26.*t153;
			t160 = t29.*t153;
			t161 = t32.*t153;
			t173 = -t165;
			t184 = t25.*t178;
			t186 = t31.*t178;
			t207 = t163+t164;
			t218 = t180+t185;
			t221 = t183+t187;
			t226 = -t24.*(t182-t183);
			t227 = -t30.*(t182-t183);
			t230 = t192+t197;
			t231 = t196+t198;
			t238 = -t5.*t6.*(t193-t196);
			t240 = t5.*t6.*(t193-t196);
			t274 = -t272;
			t277 = t133+t264;
			t280 = t108+t117+t259;
			t166 = -t158;
			t167 = -t161;
			t169 = t25.*t159;
			t170 = t25.*t161;
			t171 = t31.*t159;
			t172 = t31.*t161;
			t188 = -t186;
			t191 = t105+t159;
			t194 = t109+t161;
			t208 = t162+t173;
			t210 = t24.*t207;
			t211 = t30.*t207;
			t219 = t181+t184;
			t222 = t24.*t218;
			t223 = t30.*t218;
			t234 = t23.*t230;
			t235 = t29.*t230;
			t285 = t157+t274;
			t174 = t31.*t167;
			t195 = l22+t104+t167;
			t199 = t25.*t191;
			t200 = t31.*t191;
			t201 = -t25.*(t104+t167);
			t202 = -t31.*(t104+t167);
			t206 = t31.*(t104+t167);
			t209 = t170+t171;
			t213 = t24.*t208;
			t214 = t30.*t208;
			t215 = -t211;
			t220 = t179+t188;
			t224 = t30.*t219;
			t228 = -t223;
			t261 = t112+t235;
			t263 = t122+t234;
			t268 = -t6.*t12.*(t115-t234);
			t282 = t222+t227;
			t290 = -t23.*(t223+t24.*(t182-t183));
			t291 = -t5.*t6.*(t223+t24.*(t182-t183));
			t292 = -t29.*(t223+t24.*(t182-t183));
			t294 = t5.*t6.*(t223+t24.*(t182-t183));
			t300 = -t13.*(t113+t23.*(t223+t24.*(t182-t183)));
			t203 = -t200;
			t204 = t25.*t195;
			t205 = t31.*t195;
			t212 = t169+t174;
			t216 = t24.*t209;
			t225 = t24.*t220;
			t229 = -t224;
			t239 = t200+t201;
			t242 = t199+t206;
			t266 = t13.*t261;
			t270 = t210+t214;
			t271 = t213+t215;
			t279 = -t13.*t29.*(t211-t213);
			t281 = -t6.*t12.*t23.*(t211-t213);
			t284 = t226+t228;
			t289 = t13.*t29.*t282;
			t293 = t6.*t12.*t23.*t282;
			t297 = t118+t292;
			t298 = t120+t290;
			t217 = t30.*t212;
			t241 = t199+t205;
			t243 = t203+t204;
			t245 = t24.*t239;
			t246 = t30.*t239;
			t250 = t24.*t242;
			t251 = t30.*t242;
			t253 = -t24.*(t200-t204);
			t254 = -t30.*(t200-t204);
			t273 = t29.*t270;
			t275 = t23.*t270;
			t276 = t5.*t6.*t270;
			t283 = t225+t229;
			t299 = t6.*t12.*t297;
			t324 = t240+t266+t268;
			t332 = t289+t293+t294;
			t244 = l21+t243;
			t247 = t24.*t241;
			t248 = t30.*t241;
			t252 = -t246;
			t278 = t216+t217;
			t286 = t160+t275;
			t287 = t166+t273;
			t295 = -t6.*t12.*(t158-t273);
			t296 = t6.*t12.*(t158-t273);
			t301 = -t299;
			t302 = t245+t251;
			t309 = -t5.*t6.*(t246-t250);
			t311 = -t13.*t23.*(t246-t250);
			t312 = t5.*t6.*(t246-t250);
			t317 = -t6.*t12.*t29.*(t246-t250);
			t322 = t6.*t12.*t29.*(t246-t250);
			t327 = t276+t279+t281;
			t328 = q2w5__dt_1_.*(-t276+t13.*t29.*(t211-t213)+t6.*t12.*t23.*(t211-t213)).*-2.0;
			t333 = q2w6__dt_1_.*t332.*2.0;
			t255 = -t248;
			t256 = t24.*t244;
			t257 = t30.*t244;
			t288 = t13.*t286;
			t303 = t250+t252;
			t304 = t247+t254;
			t308 = t13.*t29.*t302;
			t313 = t6.*t12.*t23.*t302;
			t318 = -t13.*t29.*(t248+t24.*(t200-t204));
			t320 = -t6.*t12.*t23.*(t248+t24.*(t200-t204));
			t326 = t300+t301;
			t330 = t311+t322;
			t334 = -t333;
			t258 = -t256;
			t305 = t247+t257;
			t306 = t253+t255;
			t310 = t5.*t6.*t304;
			t314 = t13.*t23.*t304;
			t316 = t6.*t12.*t29.*t304;
			t325 = t288+t296;
			t335 = t308+t312+t313;
			t307 = t248+t258;
			t315 = t13.*t23.*t305;
			t319 = t6.*t12.*t29.*t305;
			t321 = -t316;
			t336 = q2w4__dt_1_.*t335.*2.0;
			t337 = t310+t318+t320;
			t338 = q2w3__dt_1_.*(-t310+t13.*t29.*(t248+t24.*(t200-t204))+t6.*t12.*t23.*(t248+t24.*(t200-t204))).*-2.0;
			t323 = -t319;
			t329 = t314+t321;
			t331 = t315+t323;
			et1 = q2w3__dt_1_.*(t333-t336+t338+q2w1__dt_1_.*t329.*2.0+q2w5__dt_1_.*(-t276+t13.*t29.*(t211-t213)+t6.*t12.*t23.*(t211-t213)).*2.0-q2w2__dt_1_.*(-t310+t13.*t29.*(t248+t24.*(t200-t204))+t6.*t12.*t23.*(t248+t24.*(t200-t204))).*2.0)+q2w6__dt_1_.*(q2w2__dt_1_.*t332.*2.0+q2w3__dt_1_.*t332.*2.0+q2w4__dt_1_.*t332.*2.0+q2w5__dt_1_.*t324.*2.0+q2w1__dt_1_.*(t299+t13.*(t113+t23.*(t223+t24.*(t182-t183)))).*2.0-q2w6__dt_1_.*(t13.*(t29.*(t224-t225)+t23.*t27.*t91)-t5.*t6.*(t24.*t219+t30.*t220)+t6.*t12.*(t23.*(t224-t225)-t27.*t29.*t91)).*2.0);
			et2 = q2w2__dt_1_.*(t333-t336+t338+q2w1__dt_1_.*t331.*2.0+q2w5__dt_1_.*(-t276+t13.*t29.*(t211-t213)+t6.*t12.*t23.*(t211-t213)).*2.0-q2w2__dt_1_.*(-t5.*t6.*t305+t13.*t29.*t307+t6.*t12.*t23.*t307).*2.0)+q1dev__dt_1_.*(q1flex__dt_1_.*t285.*-2.0+q1sup__dt_1_.*t277.*2.0+q1dev__dt_1_.*(t15.*(t50+t90-t97)-t8.*t14.*(t55-t98+l12y.*(-t73+t78+t102+t168-t74.*t78.*2.0))).*2.0);
			et3 = q2w5__dt_1_.*(q2w1__dt_1_.*t325.*-2.0+q2w6__dt_1_.*t324.*2.0+q2w2__dt_1_.*(-t276+t13.*t29.*(t211-t213)+t6.*t12.*t23.*(t211-t213)).*2.0+q2w3__dt_1_.*(-t276+t13.*t29.*(t211-t213)+t6.*t12.*t23.*(t211-t213)).*2.0+q2w4__dt_1_.*(-t276+t13.*t29.*(t211-t213)+t6.*t12.*t23.*(t211-t213)).*2.0+q2w5__dt_1_.*(-t13.*(t155+t29.*t278)+t5.*t6.*(t24.*t212-t30.*t209)+t6.*t12.*(t156-t23.*t278)).*2.0)-q1sup__dt_1_.*(q1dev__dt_1_.*t277.*-2.0+q1flex__dt_1_.*t280.*2.0+q1sup__dt_1_.*(-t15.*(t50+l12y.*t62+l12z.*t61)+t8.*t14.*(t55+l12z.*t63+l12y.*t168)+t7.*t8.*(t51+t57-l12x.*t11)).*2.0);
			et4 = q2w1__dt_1_.*(q2w2__dt_1_.*t331.*2.0+q2w3__dt_1_.*t329.*2.0-q2w5__dt_1_.*t325.*2.0+q2w6__dt_1_.*(t299+t13.*(t113+t23.*(t223+t24.*(t182-t183)))).*2.0-q2w1__dt_1_.*(t13.*(t155+t29.*t307)-t6.*t12.*(t156-t23.*t307)).*2.0-q2w4__dt_1_.*(t317+t13.*t23.*(t246-t250)).*2.0)-q2w4__dt_1_.*(t328+t334+t336+q2w2__dt_1_.*t335.*2.0+q2w3__dt_1_.*t335.*2.0+q2w1__dt_1_.*(t317+t13.*t23.*(t246-t250)).*2.0)-q1flex__dt_1_.*(q1dev__dt_1_.*t285.*2.0+q1sup__dt_1_.*t280.*2.0-q1flex__dt_1_.*(t15.*(t90-t97)-t7.*t8.*(t51+t57)+t8.*t14.*(t98-l12y.*(-t74+t79+t102+t168-t73.*t79.*2.0))).*2.0);
			out1 = et1+et2+et3+et4;
			end
			function out1 = ConsCor_4(t,in2,in3,in4,in5,in6)
			%ConsCor_4
			%    OUT1 = ConsCor_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:26
			ax__FRAME_5__dt_1_ = in6(7);
			ax__FRAME_20__dt_1_ = in6(20);
			ay__FRAME_5__dt_1_ = in6(8);
			ay__FRAME_20__dt_1_ = in6(21);
			az__FRAME_5__dt_1_ = in6(9);
			az__FRAME_20__dt_1_ = in6(22);
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			t2 = (ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t3 = (ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t4 = (ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t5 = (ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t6 = (ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t7 = (ax__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t8 = (ax__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t9 = (ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t10 = (ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t11 = (ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t12 = (ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t13 = (ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t14 = (ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t15 = (ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t16 = (ax__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t17 = (ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t18 = (ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t19 = (ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t20 = (ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t21 = (ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t22 = (ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t23 = (ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t24 = (ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t25 = (ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t26 = (ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t27 = (ay__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t28 = (ay__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t29 = (ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t30 = (ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t31 = (ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t32 = (ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t33 = (ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t34 = (az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t35 = (az__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t36 = (az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t37 = (az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t38 = (az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t39 = (az__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t40 = (az__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t41 = (az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t42 = (az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t43 = (az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t44 = (az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t45 = (az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t46 = (az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t47 = (az__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t48 = (az__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t49 = (az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t50 = -t4;
			t51 = -t6;
			t52 = -t7;
			t53 = -t8;
			t54 = -t12;
			t55 = -t13;
			t56 = -t14;
			t57 = -t17;
			t58 = -t21;
			t59 = -t24;
			t60 = -t26;
			t61 = -t31;
			t62 = -t33;
			t63 = -t35;
			t64 = -t38;
			t65 = -t41;
			t66 = -t44;
			t67 = -t46;
			t68 = t5+t25+t45;
			t69 = t10+t30+t49;
			t70 = t9+t34+t58;
			t71 = t18+t37+t55;
			t72 = t2+t29+t65;
			t73 = t14+t38+t60;
			t74 = t22+t42+t57;
			t75 = t6+t33+t67;
			t76 = (rqw__FRAME_20__dt_0_.*t68)./2.0;
			t77 = (rqy__FRAME_20__dt_0_.*t68)./2.0;
			t78 = (rqw__FRAME_5__dt_0_.*t69)./2.0;
			t79 = (rqz__FRAME_20__dt_0_.*t68)./2.0;
			t80 = (rqy__FRAME_5__dt_0_.*t69)./2.0;
			t81 = (rqz__FRAME_5__dt_0_.*t69)./2.0;
			t104 = t28+t48+t53+t68;
			t105 = t27+t47+t52+t69;
			t107 = t11+t23+t26+t56+t63+t64;
			t109 = t3+t43+t46+t51+t61+t62;
			t82 = (rqw__FRAME_20__dt_0_.*t70)./2.0;
			t83 = (rqw__FRAME_20__dt_0_.*t71)./2.0;
			t84 = (rqx__FRAME_20__dt_0_.*t70)./2.0;
			t85 = (rqx__FRAME_20__dt_0_.*t71)./2.0;
			t86 = (rqw__FRAME_5__dt_0_.*t73)./2.0;
			t87 = -t76;
			t88 = (rqx__FRAME_20__dt_0_.*t72)./2.0;
			t89 = (rqy__FRAME_20__dt_0_.*t71)./2.0;
			t90 = (rqz__FRAME_20__dt_0_.*t70)./2.0;
			t91 = (rqw__FRAME_5__dt_0_.*t74)./2.0;
			t92 = (rqx__FRAME_5__dt_0_.*t73)./2.0;
			t93 = (rqy__FRAME_20__dt_0_.*t72)./2.0;
			t94 = (rqx__FRAME_5__dt_0_.*t74)./2.0;
			t95 = (rqz__FRAME_20__dt_0_.*t72)./2.0;
			t96 = -t78;
			t97 = (rqx__FRAME_5__dt_0_.*t75)./2.0;
			t98 = (rqy__FRAME_5__dt_0_.*t74)./2.0;
			t99 = (rqz__FRAME_5__dt_0_.*t73)./2.0;
			t100 = (rqy__FRAME_5__dt_0_.*t75)./2.0;
			t101 = (rqz__FRAME_5__dt_0_.*t75)./2.0;
			t102 = t16+t20+t40+t71;
			t103 = t15+t19+t39+t74;
			t106 = t36+t54+t59+t70;
			t108 = t32+t50+t66+t72;
			et1 = -ay__FRAME_5__dt_1_.*(-t79+t81-t82-t85+t86+t93+t94-t100-(rqw__FRAME_5__dt_0_.*t107)./2.0+(rqx__FRAME_5__dt_0_.*t103)./2.0+(rqy__FRAME_5__dt_0_.*t109)./2.0+(rqz__FRAME_5__dt_0_.*t105)./2.0)+ay__FRAME_20__dt_1_.*(t79-t81+t82+t85-t86-t93-t94+t100+(rqw__FRAME_20__dt_0_.*t106)./2.0+(rqx__FRAME_20__dt_0_.*t102)./2.0-(rqy__FRAME_20__dt_0_.*t108)./2.0+(rqz__FRAME_20__dt_0_.*t104)./2.0)+az__FRAME_5__dt_1_.*(-t77+t80-t83+t84+t91-t92-t95+t101+(rqw__FRAME_5__dt_0_.*t103)./2.0+(rqx__FRAME_5__dt_0_.*t107)./2.0+(rqy__FRAME_5__dt_0_.*t105)./2.0-(rqz__FRAME_5__dt_0_.*t109)./2.0)-az__FRAME_20__dt_1_.*(t77-t80+t83-t84-t91+t92+t95-t101+(rqw__FRAME_20__dt_0_.*t102)./2.0-(rqx__FRAME_20__dt_0_.*t106)./2.0+(rqy__FRAME_20__dt_0_.*t104)./2.0+(rqz__FRAME_20__dt_0_.*t108)./2.0);
			et2 = ax__FRAME_5__dt_1_.*(t87+t88+t89+t90+t96+t97+t98+t99-(rqw__FRAME_5__dt_0_.*t105)./2.0-(rqx__FRAME_5__dt_0_.*t109)./2.0+(rqy__FRAME_5__dt_0_.*t103)./2.0-(rqz__FRAME_5__dt_0_.*t107)./2.0)-ax__FRAME_20__dt_1_.*(t87+t88+t89+t90+t96+t97+t98+t99-(rqw__FRAME_20__dt_0_.*t104)./2.0+(rqx__FRAME_20__dt_0_.*t108)./2.0+(rqy__FRAME_20__dt_0_.*t102)./2.0+(rqz__FRAME_20__dt_0_.*t106)./2.0);
			out1 = et1+et2;
			end
			function out1 = ConsCor_5(t,in2,in3,in4,in5,in6)
			%ConsCor_5
			%    OUT1 = ConsCor_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:27
			ax__FRAME_5__dt_1_ = in6(7);
			ax__FRAME_20__dt_1_ = in6(20);
			ay__FRAME_5__dt_1_ = in6(8);
			ay__FRAME_20__dt_1_ = in6(21);
			az__FRAME_5__dt_1_ = in6(9);
			az__FRAME_20__dt_1_ = in6(22);
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			t2 = (ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t3 = (ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t4 = (ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t5 = (ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t6 = (ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t7 = (ax__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t8 = (ax__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t9 = (ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t10 = (ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t11 = (ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t12 = (ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t13 = (ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t14 = (ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t15 = (ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t16 = (ax__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t17 = (ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t18 = (ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t19 = (ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t20 = (ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t21 = (ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t22 = (ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t23 = (ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t24 = (ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t25 = (ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t26 = (ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t27 = (ay__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t28 = (ay__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t29 = (ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t30 = (ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t31 = (ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t32 = (ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t33 = (ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t34 = (az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t35 = (az__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t36 = (az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t37 = (az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t38 = (az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t39 = (az__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t40 = (az__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t41 = (az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t42 = (az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t43 = (az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t44 = (az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t45 = (az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t46 = (az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t47 = (az__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t48 = (az__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t49 = (az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t50 = -t13;
			t51 = -t17;
			t52 = -t18;
			t53 = -t21;
			t54 = -t22;
			t55 = -t26;
			t56 = -t27;
			t57 = -t28;
			t58 = -t31;
			t59 = -t32;
			t60 = -t37;
			t61 = -t39;
			t62 = -t40;
			t63 = -t41;
			t64 = -t42;
			t65 = -t43;
			t66 = -t44;
			t67 = -t46;
			t68 = t5+t25+t45;
			t69 = t10+t30+t49;
			t70 = t9+t34+t53;
			t71 = t18+t37+t50;
			t72 = t2+t29+t63;
			t73 = t14+t38+t55;
			t74 = t22+t42+t51;
			t75 = t6+t33+t67;
			t76 = (rqw__FRAME_20__dt_0_.*t68)./2.0;
			t77 = (rqx__FRAME_20__dt_0_.*t68)./2.0;
			t78 = (rqw__FRAME_5__dt_0_.*t69)./2.0;
			t79 = (rqz__FRAME_20__dt_0_.*t68)./2.0;
			t80 = (rqx__FRAME_5__dt_0_.*t69)./2.0;
			t81 = (rqz__FRAME_5__dt_0_.*t69)./2.0;
			t104 = t8+t48+t57+t68;
			t105 = t7+t47+t56+t69;
			t106 = t13+t16+t20+t52+t60+t62;
			t108 = t15+t17+t19+t54+t61+t64;
			t82 = (rqw__FRAME_20__dt_0_.*t70)./2.0;
			t83 = (rqw__FRAME_20__dt_0_.*t72)./2.0;
			t84 = (rqx__FRAME_20__dt_0_.*t71)./2.0;
			t85 = (rqy__FRAME_20__dt_0_.*t70)./2.0;
			t86 = (rqw__FRAME_5__dt_0_.*t73)./2.0;
			t87 = -t76;
			t88 = (rqx__FRAME_20__dt_0_.*t72)./2.0;
			t89 = (rqy__FRAME_20__dt_0_.*t71)./2.0;
			t90 = (rqz__FRAME_20__dt_0_.*t70)./2.0;
			t91 = (rqy__FRAME_20__dt_0_.*t72)./2.0;
			t92 = (rqz__FRAME_20__dt_0_.*t71)./2.0;
			t93 = (rqw__FRAME_5__dt_0_.*t75)./2.0;
			t94 = (rqx__FRAME_5__dt_0_.*t74)./2.0;
			t95 = (rqy__FRAME_5__dt_0_.*t73)./2.0;
			t96 = -t78;
			t97 = (rqx__FRAME_5__dt_0_.*t75)./2.0;
			t98 = (rqy__FRAME_5__dt_0_.*t74)./2.0;
			t99 = (rqz__FRAME_5__dt_0_.*t73)./2.0;
			t100 = (rqy__FRAME_5__dt_0_.*t75)./2.0;
			t101 = (rqz__FRAME_5__dt_0_.*t74)./2.0;
			t102 = t12+t24+t36+t70;
			t103 = t11+t23+t35+t73;
			t107 = t4+t59+t66+t72;
			t109 = t3+t58+t65+t75;
			et1 = -ax__FRAME_5__dt_1_.*(t79-t81+t82+t84-t86-t91-t94+t100-(rqw__FRAME_5__dt_0_.*t103)./2.0+(rqx__FRAME_5__dt_0_.*t108)./2.0+(rqy__FRAME_5__dt_0_.*t109)./2.0-(rqz__FRAME_5__dt_0_.*t105)./2.0)-ax__FRAME_20__dt_1_.*(t79-t81+t82+t84-t86-t91-t94+t100+(rqw__FRAME_20__dt_0_.*t102)./2.0-(rqx__FRAME_20__dt_0_.*t106)./2.0-(rqy__FRAME_20__dt_0_.*t107)./2.0+(rqz__FRAME_20__dt_0_.*t104)./2.0)-az__FRAME_5__dt_1_.*(-t77+t80-t83-t85+t92+t93+t95-t101+(rqw__FRAME_5__dt_0_.*t109)./2.0+(rqx__FRAME_5__dt_0_.*t105)./2.0+(rqy__FRAME_5__dt_0_.*t103)./2.0+(rqz__FRAME_5__dt_0_.*t108)./2.0)+az__FRAME_20__dt_1_.*(t77-t80+t83+t85-t92-t93-t95+t101+(rqw__FRAME_20__dt_0_.*t107)./2.0+(rqx__FRAME_20__dt_0_.*t104)./2.0+(rqy__FRAME_20__dt_0_.*t102)./2.0+(rqz__FRAME_20__dt_0_.*t106)./2.0);
			et2 = ay__FRAME_5__dt_1_.*(t87+t88+t89+t90+t96+t97+t98+t99-(rqw__FRAME_5__dt_0_.*t105)./2.0+(rqx__FRAME_5__dt_0_.*t109)./2.0-(rqy__FRAME_5__dt_0_.*t108)./2.0+(rqz__FRAME_5__dt_0_.*t103)./2.0)-ay__FRAME_20__dt_1_.*(t87+t88+t89+t90+t96+t97+t98+t99-(rqw__FRAME_20__dt_0_.*t104)./2.0+(rqx__FRAME_20__dt_0_.*t107)./2.0-(rqy__FRAME_20__dt_0_.*t106)./2.0+(rqz__FRAME_20__dt_0_.*t102)./2.0);
			out1 = et1+et2;
			end
			function out1 = ConsCor_6(t,in2,in3,in4,in5,in6)
			%ConsCor_6
			%    OUT1 = ConsCor_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:28
			ax__FRAME_5__dt_1_ = in6(7);
			ax__FRAME_20__dt_1_ = in6(20);
			ay__FRAME_5__dt_1_ = in6(8);
			ay__FRAME_20__dt_1_ = in6(21);
			az__FRAME_5__dt_1_ = in6(9);
			az__FRAME_20__dt_1_ = in6(22);
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			t2 = (ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t3 = (ax__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t4 = (ax__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t5 = (ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t6 = (ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t7 = (ax__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t8 = (ax__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t9 = (ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t10 = (ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t11 = (ax__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t12 = (ax__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t13 = (ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t14 = (ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t15 = (ax__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t16 = (ax__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t17 = (ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t18 = (ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t19 = (ay__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t20 = (ay__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t21 = (ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t22 = (ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t23 = (ay__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t24 = (ay__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t25 = (ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t26 = (ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t27 = (ay__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t28 = (ay__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t29 = (ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t30 = (ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t31 = (ay__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t32 = (ay__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t33 = (ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t34 = (az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t35 = (az__FRAME_5__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t36 = (az__FRAME_20__dt_1_.*rqw__FRAME_5__dt_0_)./2.0;
			t37 = (az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t38 = (az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0;
			t39 = (az__FRAME_5__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t40 = (az__FRAME_20__dt_1_.*rqx__FRAME_5__dt_0_)./2.0;
			t41 = (az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t42 = (az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0;
			t43 = (az__FRAME_5__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t44 = (az__FRAME_20__dt_1_.*rqy__FRAME_5__dt_0_)./2.0;
			t45 = (az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t46 = (az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0;
			t47 = (az__FRAME_5__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t48 = (az__FRAME_20__dt_1_.*rqz__FRAME_5__dt_0_)./2.0;
			t49 = (az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0;
			t50 = -t13;
			t51 = -t17;
			t52 = -t18;
			t53 = -t19;
			t54 = -t20;
			t55 = -t21;
			t56 = -t22;
			t57 = -t23;
			t58 = -t24;
			t59 = -t26;
			t60 = -t35;
			t61 = -t36;
			t62 = -t37;
			t63 = -t41;
			t64 = -t42;
			t65 = -t46;
			t66 = -t47;
			t67 = -t48;
			t68 = t5+t25+t45;
			t69 = t10+t30+t49;
			t70 = t9+t34+t55;
			t71 = t18+t37+t50;
			t72 = t2+t29+t63;
			t73 = t14+t38+t59;
			t74 = t22+t42+t51;
			t75 = t6+t33+t65;
			t76 = (rqw__FRAME_20__dt_0_.*t68)./2.0;
			t77 = (rqx__FRAME_20__dt_0_.*t68)./2.0;
			t78 = (rqy__FRAME_20__dt_0_.*t68)./2.0;
			t79 = (rqw__FRAME_5__dt_0_.*t69)./2.0;
			t80 = (rqx__FRAME_5__dt_0_.*t69)./2.0;
			t81 = (rqy__FRAME_5__dt_0_.*t69)./2.0;
			t103 = t8+t28+t67+t68;
			t105 = t7+t27+t66+t69;
			t107 = t13+t16+t40+t52+t54+t62;
			t109 = t15+t17+t39+t53+t56+t64;
			t82 = (rqw__FRAME_20__dt_0_.*t71)./2.0;
			t83 = (rqx__FRAME_20__dt_0_.*t70)./2.0;
			t84 = (rqw__FRAME_20__dt_0_.*t72)./2.0;
			t85 = (rqy__FRAME_20__dt_0_.*t70)./2.0;
			t86 = -t76;
			t87 = (rqx__FRAME_20__dt_0_.*t72)./2.0;
			t88 = (rqy__FRAME_20__dt_0_.*t71)./2.0;
			t89 = (rqz__FRAME_20__dt_0_.*t70)./2.0;
			t90 = (rqw__FRAME_5__dt_0_.*t74)./2.0;
			t91 = (rqx__FRAME_5__dt_0_.*t73)./2.0;
			t92 = (rqz__FRAME_20__dt_0_.*t71)./2.0;
			t93 = (rqw__FRAME_5__dt_0_.*t75)./2.0;
			t94 = (rqy__FRAME_5__dt_0_.*t73)./2.0;
			t95 = (rqz__FRAME_20__dt_0_.*t72)./2.0;
			t96 = -t79;
			t97 = (rqx__FRAME_5__dt_0_.*t75)./2.0;
			t98 = (rqy__FRAME_5__dt_0_.*t74)./2.0;
			t99 = (rqz__FRAME_5__dt_0_.*t73)./2.0;
			t100 = (rqz__FRAME_5__dt_0_.*t74)./2.0;
			t101 = (rqz__FRAME_5__dt_0_.*t75)./2.0;
			t102 = t4+t32+t44+t72;
			t104 = t3+t31+t43+t75;
			t106 = t12+t58+t61+t70;
			t108 = t11+t57+t60+t73;
			et1 = ax__FRAME_5__dt_1_.*(t78-t81+t82-t83-t90+t91+t95-t101+(rqw__FRAME_5__dt_0_.*t109)./2.0+(rqx__FRAME_5__dt_0_.*t108)./2.0-(rqy__FRAME_5__dt_0_.*t105)./2.0-(rqz__FRAME_5__dt_0_.*t104)./2.0)+ax__FRAME_20__dt_1_.*(t78-t81+t82-t83-t90+t91+t95-t101-(rqw__FRAME_20__dt_0_.*t107)./2.0-(rqx__FRAME_20__dt_0_.*t106)./2.0+(rqy__FRAME_20__dt_0_.*t103)./2.0+(rqz__FRAME_20__dt_0_.*t102)./2.0)+ay__FRAME_5__dt_1_.*(-t77+t80-t84-t85+t92+t93+t94-t100+(rqw__FRAME_5__dt_0_.*t104)./2.0+(rqx__FRAME_5__dt_0_.*t105)./2.0+(rqy__FRAME_5__dt_0_.*t108)./2.0+(rqz__FRAME_5__dt_0_.*t109)./2.0)-ay__FRAME_20__dt_1_.*(t77-t80+t84+t85-t92-t93-t94+t100+(rqw__FRAME_20__dt_0_.*t102)./2.0+(rqx__FRAME_20__dt_0_.*t103)./2.0+(rqy__FRAME_20__dt_0_.*t106)./2.0+(rqz__FRAME_20__dt_0_.*t107)./2.0);
			et2 = az__FRAME_5__dt_1_.*(t86+t87+t88+t89+t96+t97+t98+t99-(rqw__FRAME_5__dt_0_.*t105)./2.0+(rqx__FRAME_5__dt_0_.*t104)./2.0-(rqy__FRAME_5__dt_0_.*t109)./2.0+(rqz__FRAME_5__dt_0_.*t108)./2.0)-az__FRAME_20__dt_1_.*(t86+t87+t88+t89+t96+t97+t98+t99-(rqw__FRAME_20__dt_0_.*t103)./2.0+(rqx__FRAME_20__dt_0_.*t102)./2.0-(rqy__FRAME_20__dt_0_.*t107)./2.0+(rqz__FRAME_20__dt_0_.*t106)./2.0);
			out1 = et1+et2;
			end
			function out1 = ConsCor_7(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsCor_7
			%    OUT1 = ConsCor_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:30
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
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
			%    14-Apr-2023 16:21:17
			c = in3(3,:);
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l11x = in3(7,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			l21x = in3(16,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_0_ = in2(4,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_0_ = in2(5,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_0_ = in2(6,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_0_ = in2(7,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_0_ = in2(8,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_0_ = in2(9,:);
			q2w6__dt_1_ = in2(18,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = cos(rf1Devx);
			t9 = cos(rf1Devy);
			t10 = cos(rf1Devz);
			t11 = sin(q1dev__dt_0_);
			t12 = sin(q1flex__dt_0_);
			t13 = sin(q1sup__dt_0_);
			t14 = sin(rf2Medx);
			t15 = sin(rf2Medy);
			t16 = sin(rf2Medz);
			t17 = sin(rf1Devx);
			t18 = sin(rf1Devy);
			t19 = sin(rf1Devz);
			t20 = l25+l13z;
			t21 = q2w1__dt_0_+riw1;
			t22 = q2w2__dt_0_+riw2;
			t23 = q2w3__dt_0_+riw3;
			t24 = q2w4__dt_0_+riw4;
			t25 = q2w5__dt_0_+riw5;
			t26 = q2w6__dt_0_+riw6;
			t47 = q1dev__dt_0_./2.0;
			t48 = q1flex__dt_0_./2.0;
			t49 = q1sup__dt_0_./2.0;
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = cos(t23);
			t30 = cos(t24);
			t31 = cos(t25);
			t32 = cos(t26);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = sin(t23);
			t36 = sin(t24);
			t37 = sin(t25);
			t38 = sin(t26);
			t39 = t2.*t3;
			t40 = t2.*t12;
			t41 = t3.*t11;
			t42 = t5.*t16;
			t43 = t8.*t19;
			t44 = t11.*t12;
			t45 = t14.*t16;
			t46 = t17.*t19;
			t50 = cos(t47);
			t51 = cos(t48);
			t52 = cos(t49);
			t53 = sin(t47);
			t54 = sin(t48);
			t55 = sin(t49);
			t58 = l12x.*t2.*t4;
			t61 = l12x.*t4.*t11;
			t62 = t7.*t14.*t15;
			t63 = t10.*t17.*t18;
			t66 = t5.*t7.*t15;
			t67 = t8.*t10.*t18;
			t56 = l13x.*t32;
			t57 = l13y.*t32;
			t59 = l13x.*t38;
			t60 = l13y.*t38;
			t64 = t13.*t44;
			t65 = t13.*t39;
			t68 = t13.*t40;
			t69 = t13.*t41;
			t70 = t20.*t31;
			t71 = t20.*t37;
			t73 = t50.^2;
			t74 = t51.^2;
			t75 = -t62;
			t76 = -t63;
			t82 = t45+t66;
			t83 = t46+t67;
			t72 = -t60;
			t77 = -t68;
			t78 = -t69;
			t79 = t57+t59;
			t80 = t39+t64;
			t81 = t44+t65;
			t89 = t42+t75;
			t90 = t43+t76;
			t84 = l23+t79;
			t85 = l12z.*t81;
			t86 = t56+t72;
			t87 = t40+t78;
			t88 = t41+t77;
			t101 = t30.*t37.*t79;
			t102 = t36.*t37.*t79;
			t91 = l24+t86;
			t92 = l12z.*t87;
			t94 = t30.*t84;
			t95 = t36.*t84;
			t96 = t30.*t86;
			t97 = t36.*t86;
			t104 = -t102;
			t93 = -t92;
			t98 = -t94;
			t99 = t31.*t91;
			t100 = t37.*t91;
			t118 = t97+t101;
			t119 = t96+t104;
			t103 = -t100;
			t105 = t71+t99;
			t120 = t29.*t118;
			t121 = t35.*t118;
			t122 = t29.*t119;
			t123 = t35.*t119;
			t106 = t70+t103;
			t107 = t27.*t105;
			t108 = t33.*t105;
			t112 = t29.*t30.*t105;
			t113 = t29.*t36.*t105;
			t114 = t30.*t35.*t105;
			t115 = t35.*t36.*t105;
			t124 = -t121;
			t140 = t120+t123;
			t143 = -t28.*(t121-t122);
			t109 = t30.*t106;
			t110 = -t108;
			t111 = t36.*t106;
			t117 = -t115;
			t136 = t113+t114;
			t141 = t122+t124;
			t142 = t34.*t140;
			t116 = -t111;
			t125 = t95+t109;
			t126 = t98+t111;
			t137 = t112+t117;
			t138 = t28.*t136;
			t144 = -t142;
			t127 = l22+t94+t116;
			t128 = t29.*t125;
			t129 = t35.*t125;
			t130 = -t29.*(t94+t116);
			t131 = -t35.*(t94+t116);
			t135 = t35.*(t94+t116);
			t139 = t34.*t137;
			t160 = t143+t144;
			t132 = -t129;
			t133 = t29.*t127;
			t134 = t35.*t127;
			t145 = t129+t130;
			t147 = t128+t135;
			t159 = t138+t139;
			t146 = t128+t134;
			t148 = t132+t133;
			t150 = t34.*t145;
			t153 = t28.*t147;
			t155 = -t34.*(t129-t133);
			t149 = l21+t148;
			t151 = t28.*t146;
			t152 = t34.*t146;
			t154 = -t150;
			t156 = t28.*t149;
			t157 = t34.*t149;
			t161 = t153+t154;
			t162 = t151+t155;
			t158 = -t156;
			t163 = t151+t157;
			t164 = t152+t158;
			t165 = t27.*t164;
			t166 = t33.*t164;
			t167 = t107+t166;
			t168 = t110+t165;
			et1 = -q1flex__dt_1_.*(t90.*(l12z.*t80-l12y.*(t51.*t54.*2.0-t51.*t54.*t73.*4.0-t50.*t52.*t53.*t54.^2.*t55.*4.0+t50.*t52.*t53.*t55.*t74.*4.0))+t83.*(l12y.*t3.*t4-l12z.*t4.*t12)+t9.*t10.*(l12y.*t81+l12z.*t88))+q2w3__dt_1_.*(-t82.*(t152+t28.*(t129-t133))+t27.*t89.*t162+t6.*t7.*t33.*t162)-q2w4__dt_1_.*(t82.*(t28.*t145+t34.*t147)+t27.*t89.*(t150-t153)+t6.*t7.*t33.*(t150-t153))+q2w2__dt_1_.*(-t82.*t164+t27.*t89.*t163+t6.*t7.*t33.*t163)-q2w1__dt_1_.*(t89.*t167+t6.*t7.*(t108-t165));
			et2 = q1sup__dt_1_.*(t90.*(l12y.*(t50.*t51.*t52.^2.*t53.*t54.*4.0-t50.*t51.*t53.*t54.*t55.^2.*4.0)-l12x.*t11.*t13+l12z.*t4.*t41)+t83.*(l12x.*t4+l12z.*t3.*t13+l12y.*t12.*t13)-t9.*t10.*(-l12x.*t2.*t13+l12y.*t4.*t40+l12z.*t4.*t39))+q2w6__dt_1_.*(t89.*(t27.*(t142+t28.*(t121-t122))+t31.*t33.*t79)+t82.*(t28.*t140-t34.*(t121-t122))+t6.*t7.*(t33.*(t142+t28.*(t121-t122))-t27.*t31.*t79))-q2w5__dt_1_.*(t89.*(t33.*t106+t27.*t159)+t82.*(t28.*t137-t34.*t136)-t6.*t7.*(t27.*t106-t33.*t159));
			et3 = q1dev__dt_1_.*(t90.*(t58+t85+l12y.*(t50.*t53.*2.0-t50.*t53.*t74.*4.0-t51.*t52.*t53.^2.*t54.*t55.*4.0+t51.*t52.*t54.*t55.*t73.*4.0))+t9.*t10.*(t61+t93+l12y.*t80));
			out1 = -c.*(et1+et2+et3)+c.*(l11x-l21x-t82.*t163-t90.*(t61+t93+l12y.*(t73.*-2.0-t74.*2.0+t73.*t74.*4.0+t50.*t51.*t52.*t53.*t54.*t55.*8.0+1.0))+t89.*(t108-t165)+t83.*(-l12x.*t13+l12z.*t3.*t4+l12y.*t4.*t12)-t6.*t7.*t167+t9.*t10.*(t58+t85-l12y.*t88));
			end
			function out1 = ConsGF_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsGF_2
			%    OUT1 = ConsGF_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:21
			c = in3(3,:);
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l11y = in3(8,:);
			l12x = in3(10,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			l21y = in3(17,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_0_ = in2(4,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_0_ = in2(5,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_0_ = in2(6,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_0_ = in2(7,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_0_ = in2(8,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_0_ = in2(9,:);
			q2w6__dt_1_ = in2(18,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf2Medz = in3(35,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			rf1Devz = in3(29,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf2Medz);
			t8 = cos(rf1Devx);
			t9 = cos(rf1Devy);
			t10 = cos(rf1Devz);
			t11 = sin(q1dev__dt_0_);
			t12 = sin(q1flex__dt_0_);
			t13 = sin(q1sup__dt_0_);
			t14 = sin(rf2Medx);
			t15 = sin(rf2Medy);
			t16 = sin(rf2Medz);
			t17 = sin(rf1Devx);
			t18 = sin(rf1Devy);
			t19 = sin(rf1Devz);
			t20 = l25+l13z;
			t21 = q2w1__dt_0_+riw1;
			t22 = q2w2__dt_0_+riw2;
			t23 = q2w3__dt_0_+riw3;
			t24 = q2w4__dt_0_+riw4;
			t25 = q2w5__dt_0_+riw5;
			t26 = q2w6__dt_0_+riw6;
			t45 = q1dev__dt_0_./2.0;
			t46 = q1flex__dt_0_./2.0;
			t47 = q1sup__dt_0_./2.0;
			t48 = rf2Medx./2.0;
			t49 = rf2Medy./2.0;
			t50 = rf2Medz./2.0;
			t51 = rf1Devx./2.0;
			t52 = rf1Devy./2.0;
			t53 = rf1Devz./2.0;
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = cos(t23);
			t30 = cos(t24);
			t31 = cos(t25);
			t32 = cos(t26);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = sin(t23);
			t36 = sin(t24);
			t37 = sin(t25);
			t38 = sin(t26);
			t39 = t2.*t3;
			t40 = t2.*t12;
			t41 = t3.*t11;
			t42 = t7.*t14;
			t43 = t10.*t17;
			t44 = t11.*t12;
			t54 = cos(t45);
			t55 = cos(t46);
			t56 = cos(t47);
			t57 = cos(t48);
			t58 = cos(t49);
			t59 = cos(t50);
			t60 = cos(t51);
			t61 = cos(t52);
			t62 = cos(t53);
			t63 = sin(t45);
			t64 = sin(t46);
			t65 = sin(t47);
			t66 = sin(t48);
			t67 = sin(t49);
			t68 = sin(t50);
			t69 = sin(t51);
			t70 = sin(t52);
			t71 = sin(t53);
			t74 = l12x.*t2.*t4;
			t77 = l12x.*t4.*t11;
			t78 = t5.*t15.*t16;
			t79 = t8.*t18.*t19;
			t72 = l13x.*t32;
			t73 = l13y.*t32;
			t75 = l13x.*t38;
			t76 = l13y.*t38;
			t80 = t13.*t44;
			t81 = t13.*t39;
			t82 = t13.*t40;
			t83 = t13.*t41;
			t84 = t20.*t31;
			t85 = t20.*t37;
			t87 = t54.^2;
			t88 = t55.^2;
			t89 = t57.^2;
			t90 = t59.^2;
			t91 = t60.^2;
			t92 = t62.^2;
			t93 = -t78;
			t94 = -t79;
			t143 = t57.*t58.*t59.*t66.*t67.*t68.*8.0;
			t144 = t60.*t61.*t62.*t69.*t70.*t71.*8.0;
			t86 = -t76;
			t95 = -t82;
			t96 = -t83;
			t97 = t89.*2.0;
			t98 = t90.*2.0;
			t99 = t91.*2.0;
			t100 = t92.*2.0;
			t105 = t73+t75;
			t106 = t39+t80;
			t107 = t44+t81;
			t113 = t42+t93;
			t114 = t43+t94;
			t117 = t89.*t90.*4.0;
			t118 = t91.*t92.*4.0;
			t101 = -t97;
			t102 = -t98;
			t103 = -t99;
			t104 = -t100;
			t108 = l23+t105;
			t109 = l12z.*t107;
			t110 = t72+t86;
			t111 = t40+t96;
			t112 = t41+t95;
			t127 = t30.*t37.*t105;
			t128 = t36.*t37.*t105;
			t115 = l24+t110;
			t116 = l12z.*t111;
			t120 = t30.*t108;
			t121 = t36.*t108;
			t122 = t30.*t110;
			t123 = t36.*t110;
			t130 = -t128;
			t165 = t101+t102+t117+t143+1.0;
			t166 = t103+t104+t118+t144+1.0;
			t119 = -t116;
			t124 = -t120;
			t125 = t31.*t115;
			t126 = t37.*t115;
			t146 = t123+t127;
			t147 = t122+t130;
			t129 = -t126;
			t131 = t85+t125;
			t148 = t29.*t146;
			t149 = t35.*t146;
			t150 = t29.*t147;
			t151 = t35.*t147;
			t132 = t84+t129;
			t133 = t27.*t131;
			t134 = t33.*t131;
			t138 = t29.*t30.*t131;
			t139 = t29.*t36.*t131;
			t140 = t30.*t35.*t131;
			t141 = t35.*t36.*t131;
			t152 = -t149;
			t170 = t148+t151;
			t173 = -t28.*(t149-t150);
			t135 = t30.*t132;
			t136 = -t134;
			t137 = t36.*t132;
			t145 = -t141;
			t164 = t139+t140;
			t171 = t150+t152;
			t172 = t34.*t170;
			t142 = -t137;
			t153 = t121+t135;
			t154 = t124+t137;
			t167 = t138+t145;
			t168 = t28.*t164;
			t174 = -t172;
			t155 = l22+t120+t142;
			t156 = t29.*t153;
			t157 = t35.*t153;
			t158 = -t29.*(t120+t142);
			t159 = -t35.*(t120+t142);
			t163 = t35.*(t120+t142);
			t169 = t34.*t167;
			t190 = t173+t174;
			t160 = -t157;
			t161 = t29.*t155;
			t162 = t35.*t155;
			t175 = t157+t158;
			t177 = t156+t163;
			t189 = t168+t169;
			t176 = t156+t162;
			t178 = t160+t161;
			t180 = t34.*t175;
			t183 = t28.*t177;
			t185 = -t34.*(t157-t161);
			t179 = l21+t178;
			t181 = t28.*t176;
			t182 = t34.*t176;
			t184 = -t180;
			t186 = t28.*t179;
			t187 = t34.*t179;
			t191 = t183+t184;
			t192 = t181+t185;
			t188 = -t186;
			t193 = t181+t187;
			t194 = t182+t188;
			t195 = t27.*t194;
			t196 = t33.*t194;
			t197 = t133+t196;
			t198 = t136+t195;
			et1 = q1flex__dt_1_.*(t166.*(l12z.*t106-l12y.*(t55.*t64.*2.0-t55.*t64.*t87.*4.0-t54.*t56.*t63.*t64.^2.*t65.*4.0+t54.*t56.*t63.*t65.*t88.*4.0))+t114.*(l12y.*t3.*t4-l12z.*t4.*t12)-t9.*t19.*(l12y.*t107+l12z.*t112))+q2w3__dt_1_.*(t113.*(t182+t28.*(t157-t161))-t27.*t165.*t192+t6.*t16.*t33.*t192)+q2w4__dt_1_.*(t113.*(t28.*t175+t34.*t177)+t27.*t165.*(t180-t183)-t6.*t16.*t33.*(t180-t183))+q2w2__dt_1_.*(t113.*t194-t27.*t165.*t193+t6.*t16.*t33.*t193)+q2w1__dt_1_.*(t165.*t197-t6.*t16.*(t134-t195));
			et2 = -q1sup__dt_1_.*(t166.*(l12y.*(t54.*t55.*t56.^2.*t63.*t64.*4.0-t54.*t55.*t63.*t64.*t65.^2.*4.0)-l12x.*t11.*t13+l12z.*t4.*t41)+t114.*(l12x.*t4+l12z.*t3.*t13+l12y.*t12.*t13)+t9.*t19.*(-l12x.*t2.*t13+l12y.*t4.*t40+l12z.*t4.*t39))-q2w6__dt_1_.*(t165.*(t27.*(t172+t28.*(t149-t150))+t31.*t33.*t105)+t113.*(t28.*t170-t34.*(t149-t150))-t6.*t16.*(t33.*(t172+t28.*(t149-t150))-t27.*t31.*t105))+q2w5__dt_1_.*(t113.*(t28.*t167-t34.*t164)+t165.*(t33.*t132+t27.*t189)+t6.*t16.*(t27.*t132-t33.*t189));
			et3 = -q1dev__dt_1_.*(t166.*(t74+t109+l12y.*(t54.*t63.*2.0-t54.*t63.*t88.*4.0-t55.*t56.*t63.^2.*t64.*t65.*4.0+t55.*t56.*t64.*t65.*t87.*4.0))-t9.*t19.*(t77+t119+l12y.*t106));
			out1 = -c.*(et1+et2+et3)+c.*(l11y-l21y+t113.*t193+t166.*(t77+t119+l12y.*(t87.*-2.0-t88.*2.0+t87.*t88.*4.0+t54.*t55.*t56.*t63.*t64.*t65.*8.0+1.0))-t165.*(t134-t195)-t114.*(-l12x.*t13+l12z.*t3.*t4+l12y.*t4.*t12)-t6.*t16.*t197+t9.*t19.*(t74+t109-l12y.*t112));
			end
			function out1 = ConsGF_3(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsGF_3
			%    OUT1 = ConsGF_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:25
			c = in3(3,:);
			l21 = in3(19,:);
			l22 = in3(20,:);
			l23 = in3(21,:);
			l24 = in3(22,:);
			l25 = in3(23,:);
			l12x = in3(10,:);
			l11z = in3(9,:);
			l12y = in3(11,:);
			l13x = in3(13,:);
			l12z = in3(12,:);
			l13y = in3(14,:);
			l13z = in3(15,:);
			l21z = in3(18,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
			q2w1__dt_0_ = in2(4,:);
			q2w1__dt_1_ = in2(13,:);
			q2w2__dt_0_ = in2(5,:);
			q2w2__dt_1_ = in2(14,:);
			q2w3__dt_0_ = in2(6,:);
			q2w3__dt_1_ = in2(15,:);
			q2w4__dt_0_ = in2(7,:);
			q2w4__dt_1_ = in2(16,:);
			q2w5__dt_0_ = in2(8,:);
			q2w5__dt_1_ = in2(17,:);
			q2w6__dt_0_ = in2(9,:);
			q2w6__dt_1_ = in2(18,:);
			rf2Medx = in3(33,:);
			rf2Medy = in3(34,:);
			rf1Devx = in3(27,:);
			rf1Devy = in3(28,:);
			riw1 = in3(36,:);
			riw2 = in3(37,:);
			riw3 = in3(38,:);
			riw4 = in3(39,:);
			riw5 = in3(40,:);
			riw6 = in3(41,:);
			t2 = cos(q1dev__dt_0_);
			t3 = cos(q1flex__dt_0_);
			t4 = cos(q1sup__dt_0_);
			t5 = cos(rf2Medx);
			t6 = cos(rf2Medy);
			t7 = cos(rf1Devx);
			t8 = cos(rf1Devy);
			t9 = sin(q1dev__dt_0_);
			t10 = sin(q1flex__dt_0_);
			t11 = sin(q1sup__dt_0_);
			t12 = sin(rf2Medx);
			t13 = sin(rf2Medy);
			t14 = sin(rf1Devx);
			t15 = sin(rf1Devy);
			t16 = l25+l13z;
			t17 = q2w1__dt_0_+riw1;
			t18 = q2w2__dt_0_+riw2;
			t19 = q2w3__dt_0_+riw3;
			t20 = q2w4__dt_0_+riw4;
			t21 = q2w5__dt_0_+riw5;
			t22 = q2w6__dt_0_+riw6;
			t39 = q1dev__dt_0_./2.0;
			t40 = q1flex__dt_0_./2.0;
			t41 = q1sup__dt_0_./2.0;
			t23 = cos(t17);
			t24 = cos(t18);
			t25 = cos(t19);
			t26 = cos(t20);
			t27 = cos(t21);
			t28 = cos(t22);
			t29 = sin(t17);
			t30 = sin(t18);
			t31 = sin(t19);
			t32 = sin(t20);
			t33 = sin(t21);
			t34 = sin(t22);
			t35 = t2.*t3;
			t36 = t2.*t10;
			t37 = t3.*t9;
			t38 = t9.*t10;
			t42 = cos(t39);
			t43 = cos(t40);
			t44 = cos(t41);
			t45 = sin(t39);
			t46 = sin(t40);
			t47 = sin(t41);
			t50 = l12x.*t2.*t4;
			t53 = l12x.*t4.*t9;
			t48 = l13x.*t28;
			t49 = l13y.*t28;
			t51 = l13x.*t34;
			t52 = l13y.*t34;
			t54 = t11.*t38;
			t55 = t11.*t35;
			t56 = t11.*t36;
			t57 = t11.*t37;
			t58 = t16.*t27;
			t59 = t16.*t33;
			t61 = t42.^2;
			t62 = t43.^2;
			t60 = -t52;
			t63 = -t56;
			t64 = -t57;
			t65 = t49+t51;
			t66 = t35+t54;
			t67 = t38+t55;
			t68 = l23+t65;
			t69 = l12z.*t67;
			t70 = t48+t60;
			t71 = t36+t64;
			t72 = t37+t63;
			t83 = t26.*t33.*t65;
			t84 = t32.*t33.*t65;
			t73 = l24+t70;
			t74 = l12z.*t71;
			t76 = t26.*t68;
			t77 = t32.*t68;
			t78 = t26.*t70;
			t79 = t32.*t70;
			t86 = -t84;
			t75 = -t74;
			t80 = -t76;
			t81 = t27.*t73;
			t82 = t33.*t73;
			t100 = t79+t83;
			t101 = t78+t86;
			t85 = -t82;
			t87 = t59+t81;
			t102 = t25.*t100;
			t103 = t31.*t100;
			t104 = t25.*t101;
			t105 = t31.*t101;
			t88 = t58+t85;
			t89 = t23.*t87;
			t90 = t29.*t87;
			t94 = t25.*t26.*t87;
			t95 = t25.*t32.*t87;
			t96 = t26.*t31.*t87;
			t97 = t31.*t32.*t87;
			t106 = -t103;
			t122 = t102+t105;
			t125 = -t24.*(t103-t104);
			t91 = t26.*t88;
			t92 = -t90;
			t93 = t32.*t88;
			t99 = -t97;
			t118 = t95+t96;
			t123 = t104+t106;
			t124 = t30.*t122;
			t98 = -t93;
			t107 = t77+t91;
			t108 = t80+t93;
			t119 = t94+t99;
			t120 = t24.*t118;
			t126 = -t124;
			t109 = l22+t76+t98;
			t110 = t25.*t107;
			t111 = t31.*t107;
			t112 = -t25.*(t76+t98);
			t113 = -t31.*(t76+t98);
			t117 = t31.*(t76+t98);
			t121 = t30.*t119;
			t142 = t125+t126;
			t114 = -t111;
			t115 = t25.*t109;
			t116 = t31.*t109;
			t127 = t111+t112;
			t129 = t110+t117;
			t141 = t120+t121;
			t128 = t110+t116;
			t130 = t114+t115;
			t132 = t30.*t127;
			t135 = t24.*t129;
			t137 = -t30.*(t111-t115);
			t131 = l21+t130;
			t133 = t24.*t128;
			t134 = t30.*t128;
			t136 = -t132;
			t138 = t24.*t131;
			t139 = t30.*t131;
			t143 = t135+t136;
			t144 = t133+t137;
			t140 = -t138;
			t145 = t133+t139;
			t146 = t134+t140;
			t147 = t23.*t146;
			t148 = t29.*t146;
			t149 = t89+t148;
			t150 = t92+t147;
			et1 = -q1sup__dt_1_.*(t15.*(-l12x.*t2.*t11+l12y.*t4.*t36+l12z.*t4.*t35)-t8.*t14.*(l12y.*(t42.*t43.*t44.^2.*t45.*t46.*4.0-t42.*t43.*t45.*t46.*t47.^2.*4.0)-l12x.*t9.*t11+l12z.*t4.*t37)+t7.*t8.*(l12x.*t4+l12z.*t3.*t11+l12y.*t10.*t11))+q2w6__dt_1_.*(t13.*(t29.*(t124+t24.*(t103-t104))-t23.*t27.*t65)-t5.*t6.*(t24.*t122-t30.*(t103-t104))+t6.*t12.*(t23.*(t124+t24.*(t103-t104))+t27.*t29.*t65))+q2w5__dt_1_.*(t13.*(t23.*t88-t29.*t141)-t6.*t12.*(t29.*t88+t23.*t141)+t5.*t6.*(t24.*t119-t30.*t118))-q2w1__dt_1_.*(t13.*(t90-t147)+t6.*t12.*t149);
			et2 = -q1flex__dt_1_.*(t15.*(l12y.*t67+l12z.*t72)-t7.*t8.*(l12y.*t3.*t4-l12z.*t4.*t10)+t8.*t14.*(l12z.*t66-l12y.*(t43.*t46.*2.0-t43.*t46.*t61.*4.0-t42.*t44.*t45.*t46.^2.*t47.*4.0+t42.*t44.*t45.*t47.*t62.*4.0)))+q2w3__dt_1_.*(t5.*t6.*(t134+t24.*(t111-t115))+t13.*t29.*t144+t6.*t12.*t23.*t144)-q2w4__dt_1_.*(-t5.*t6.*(t24.*t127+t30.*t129)+t13.*t29.*(t132-t135)+t6.*t12.*t23.*(t132-t135))+q2w2__dt_1_.*(t5.*t6.*t146+t13.*t29.*t145+t6.*t12.*t23.*t145)+q1dev__dt_1_.*(t15.*(t53+t75+l12y.*t66)+t8.*t14.*(t50+t69+l12y.*(t42.*t45.*2.0-t42.*t45.*t62.*4.0-t43.*t44.*t45.^2.*t46.*t47.*4.0+t43.*t44.*t46.*t47.*t61.*4.0)));
			out1 = c.*(l11z-l21z+t13.*t149-t15.*(t50+t69-l12y.*t72)+t7.*t8.*(-l12x.*t11+l12z.*t3.*t4+l12y.*t4.*t10)-t5.*t6.*t145+t8.*t14.*(t53+t75+l12y.*(t61.*-2.0-t62.*2.0+t61.*t62.*4.0+t42.*t43.*t44.*t45.*t46.*t47.*8.0+1.0))-t6.*t12.*(t90-t147))+c.*(et1+et2);
			end
			function out1 = ConsGF_4(t,in2,in3,in4,in5,in6)
			%ConsGF_4
			%    OUT1 = ConsGF_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:27
			ax__FRAME_5__dt_1_ = in6(7);
			ax__FRAME_20__dt_1_ = in6(20);
			ay__FRAME_5__dt_1_ = in6(8);
			ay__FRAME_20__dt_1_ = in6(21);
			az__FRAME_5__dt_1_ = in6(9);
			az__FRAME_20__dt_1_ = in6(22);
			c = in3(3,:);
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			et1 = -c.*(rqw__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_-rqw__FRAME_20__dt_0_.*rqx__FRAME_5__dt_0_+rqy__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_-rqy__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_);
			et2 = c.*(rqw__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0)+rqy__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0)+rqx__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0)-rqz__FRAME_20__dt_0_.*(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_.*(-1.0./2.0)+(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0)-rqw__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0)-rqy__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0)-rqx__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0)+rqz__FRAME_5__dt_0_.*(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_.*(-1.0./2.0)+(ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0));
			out1 = et1+et2;
			end
			function out1 = ConsGF_5(t,in2,in3,in4,in5,in6)
			%ConsGF_5
			%    OUT1 = ConsGF_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:28
			ax__FRAME_5__dt_1_ = in6(7);
			ax__FRAME_20__dt_1_ = in6(20);
			ay__FRAME_5__dt_1_ = in6(8);
			ay__FRAME_20__dt_1_ = in6(21);
			az__FRAME_5__dt_1_ = in6(9);
			az__FRAME_20__dt_1_ = in6(22);
			c = in3(3,:);
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			et1 = -c.*(rqw__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_-rqw__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_-rqx__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_+rqx__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_);
			et2 = c.*(rqw__FRAME_20__dt_0_.*(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_.*(-1.0./2.0)+(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0)-rqx__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0)-rqw__FRAME_5__dt_0_.*(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_.*(-1.0./2.0)+(ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0)+rqx__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0)+rqy__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0)+rqz__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0)-rqy__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0)-rqz__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0));
			out1 = et1+et2;
			end
			function out1 = ConsGF_6(t,in2,in3,in4,in5,in6)
			%ConsGF_6
			%    OUT1 = ConsGF_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:29
			ax__FRAME_5__dt_1_ = in6(7);
			ax__FRAME_20__dt_1_ = in6(20);
			ay__FRAME_5__dt_1_ = in6(8);
			ay__FRAME_20__dt_1_ = in6(21);
			az__FRAME_5__dt_1_ = in6(9);
			az__FRAME_20__dt_1_ = in6(22);
			c = in3(3,:);
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			et1 = -c.*(rqw__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_-rqw__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_+rqx__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_-rqx__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_);
			et2 = c.*(rqw__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0-(ay__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0)+rqx__FRAME_20__dt_0_.*(ax__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_.*(-1.0./2.0)+(ay__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0)-rqw__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0-(ay__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0)-rqy__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqw__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0-(az__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0)-rqx__FRAME_5__dt_0_.*(ax__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_.*(-1.0./2.0)+(ay__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0)+rqz__FRAME_20__dt_0_.*((ax__FRAME_5__dt_1_.*rqx__FRAME_5__dt_0_)./2.0+(ay__FRAME_5__dt_1_.*rqy__FRAME_5__dt_0_)./2.0+(az__FRAME_5__dt_1_.*rqz__FRAME_5__dt_0_)./2.0)+rqy__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqw__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0-(az__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0)-rqz__FRAME_5__dt_0_.*((ax__FRAME_20__dt_1_.*rqx__FRAME_20__dt_0_)./2.0+(ay__FRAME_20__dt_1_.*rqy__FRAME_20__dt_0_)./2.0+(az__FRAME_20__dt_1_.*rqz__FRAME_20__dt_0_)./2.0));
			out1 = et1+et2;
			end
			function out1 = ConsGF_7(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsGF_7
			%    OUT1 = ConsGF_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:31
			c = in3(3,:);
			q1dev__dt_0_ = in2(1,:);
			q1dev__dt_1_ = in2(10,:);
			q1flex__dt_0_ = in2(2,:);
			q1flex__dt_1_ = in2(11,:);
			q1sup__dt_0_ = in2(3,:);
			q1sup__dt_1_ = in2(12,:);
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
			%    14-Apr-2023 16:21:17
			out1 = 0.0;
			end
			function out1 = ConsFVJac_2(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsFVJac_2
			%    OUT1 = ConsFVJac_2(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:22
			out1 = 0.0;
			end
			function out1 = ConsFVJac_3(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsFVJac_3
			%    OUT1 = ConsFVJac_3(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:26
			out1 = 0.0;
			end
			function out1 = ConsFVJac_4(t,in2,in3,in4,in5,in6)
			%ConsFVJac_4
			%    OUT1 = ConsFVJac_4(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:27
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			t2 = (rqw__FRAME_5__dt_0_.*rqw__FRAME_20__dt_0_)./2.0;
			t3 = (rqw__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t4 = (rqw__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
			t5 = (rqx__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
			t6 = (rqw__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t7 = (rqw__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
			t8 = (rqx__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t9 = (rqx__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
			t10 = (rqx__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t11 = (rqx__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
			t12 = (rqy__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t13 = (rqz__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t14 = -t4;
			t15 = -t6;
			t16 = -t8;
			t17 = -t10;
			t18 = t3+t11+t14+t17;
			t19 = t7+t9+t15+t16;
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t2+t5+t12+t13,-t2-t5-t12-t13,t19,t19,t18,t18],[2,6]);
			end
			function out1 = ConsFVJac_5(t,in2,in3,in4,in5,in6)
			%ConsFVJac_5
			%    OUT1 = ConsFVJac_5(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:28
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			t2 = (rqw__FRAME_5__dt_0_.*rqw__FRAME_20__dt_0_)./2.0;
			t3 = (rqw__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
			t4 = (rqw__FRAME_20__dt_0_.*rqx__FRAME_5__dt_0_)./2.0;
			t5 = (rqx__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
			t6 = (rqw__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t7 = (rqw__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
			t8 = (rqx__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t9 = (rqx__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
			t10 = (rqy__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t11 = (rqy__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t12 = (rqy__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
			t13 = (rqz__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t14 = -t3;
			t15 = -t7;
			t16 = -t9;
			t17 = -t11;
			t18 = t4+t12+t14+t17;
			t19 = t6+t8+t15+t16;
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t19,t19,t2+t5+t10+t13,-t2-t5-t10-t13,t18,t18],[2,6]);
			end
			function out1 = ConsFVJac_6(t,in2,in3,in4,in5,in6)
			%ConsFVJac_6
			%    OUT1 = ConsFVJac_6(T,IN2,IN3,IN4,IN5,IN6)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:29
			rqw__FRAME_5__dt_0_ = in6(10);
			rqw__FRAME_20__dt_0_ = in6(23);
			rqx__FRAME_5__dt_0_ = in6(11);
			rqx__FRAME_20__dt_0_ = in6(24);
			rqy__FRAME_5__dt_0_ = in6(12);
			rqy__FRAME_20__dt_0_ = in6(25);
			rqz__FRAME_5__dt_0_ = in6(13);
			rqz__FRAME_20__dt_0_ = in6(26);
			t2 = (rqw__FRAME_5__dt_0_.*rqw__FRAME_20__dt_0_)./2.0;
			t3 = (rqw__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
			t4 = (rqw__FRAME_20__dt_0_.*rqx__FRAME_5__dt_0_)./2.0;
			t5 = (rqw__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t6 = (rqw__FRAME_20__dt_0_.*rqy__FRAME_5__dt_0_)./2.0;
			t7 = (rqx__FRAME_5__dt_0_.*rqx__FRAME_20__dt_0_)./2.0;
			t8 = (rqx__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t9 = (rqx__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
			t10 = (rqy__FRAME_5__dt_0_.*rqy__FRAME_20__dt_0_)./2.0;
			t11 = (rqy__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t12 = (rqy__FRAME_20__dt_0_.*rqz__FRAME_5__dt_0_)./2.0;
			t13 = (rqz__FRAME_5__dt_0_.*rqz__FRAME_20__dt_0_)./2.0;
			t14 = -t4;
			t15 = -t5;
			t16 = -t9;
			t17 = -t12;
			t18 = t3+t11+t14+t17;
			t19 = t6+t8+t15+t16;
			out1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t19,t19,t18,t18,t2+t7+t10+t13,-t2-t7-t10-t13],[2,6]);
			end
			function out1 = ConsFVJac_7(t,in2,in3,in4,in5,SUBSVECTOR__)
			%ConsFVJac_7
			%    OUT1 = ConsFVJac_7(T,IN2,IN3,IN4,IN5,SUBSVECTOR__)
			%    This function was generated by the Symbolic Math Toolbox version 9.1.
			%    14-Apr-2023 16:21:31
			out1 = 0.0;
			end
		end
	end

end

