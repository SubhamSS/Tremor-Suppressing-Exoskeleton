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

