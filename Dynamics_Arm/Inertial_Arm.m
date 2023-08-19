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

