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
