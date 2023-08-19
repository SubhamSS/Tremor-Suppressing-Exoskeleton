function [frameTF,frameVel,frameCorAcc,frameJacobian]=frameKinematics_TAWE(frame,t,q,p,u,s,TF,Vel,CorAcc,Jacobian)
%% Author: Jiamin Wang; Updated: 2021-12-15;
    
    DOF=numel(q)/2;
    
    FramePath=0;
    thisTF=eye(4,4);
    thisTransVel=zeros(3,1);
    thisCorTransAcc=zeros(3,1);
    thisTransJacobian=zeros(3,DOF);
    thisAngVel=zeros(3,1);
    thisCorAngAcc=zeros(3,1);
    thisAngJacobian=zeros(3,DOF);
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

    for jj=numel(FramePath):-1:1
        curNum=FramePath(jj);
        curTF=TF(:,:,curNum);
        curTransVel=Vel(1:3,curNum);
        curAngVel=Vel(4:6,curNum);
        curCorTransAcc=CorAcc(1:3,curNum);
        curCorAngAcc=CorAcc(4:6,curNum);
        curTransJacobian=Jacobian(1:3,:,curNum);
        curAngJacobian=Jacobian(4:6,:,curNum);
        curRotMat=curTF(1:3,1:3);
        thisTransDis=curRotMat*thisTF(1:3,4);
        curRadVel=cross(curAngVel,thisTransDis);
        curRotAngVel=curRotMat*thisAngVel;
        thisRotTransVel=curRotMat*thisTransVel;
        thisCorTransAcc=curCorTransAcc...
                        +cross(curCorAngAcc,thisTransDis)...
                        +curRotMat*thisCorTransAcc...
                        +2*cross(curAngVel,thisRotTransVel)...
                        +cross(curAngVel,curRadVel);
        thisCorAngAcc=  curCorAngAcc...
                        +curRotMat*thisCorAngAcc...
                        +cross(curAngVel,curRotAngVel);
        thisTransVel=   curTransVel...
                        +thisRotTransVel...
                        +curRadVel;
        thisAngVel=curAngVel+curRotAngVel;
        thisTransJacobian=  curTransJacobian...
                            +curRotMat*thisTransJacobian...
                            -skew3(thisTransDis)*curAngJacobian;
        thisAngJacobian=curAngJacobian+curRotMat*thisAngJacobian;
        thisTF=curTF*thisTF;
    end
    frameTF=thisTF;
    frameVel=[thisTransVel;thisAngVel];
    frameCorAcc=[thisCorTransAcc;thisCorAngAcc];
    frameJacobian=[thisTransJacobian;thisAngJacobian];
    
    function output=skew3(w)
        output=[0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ];
    end

end

