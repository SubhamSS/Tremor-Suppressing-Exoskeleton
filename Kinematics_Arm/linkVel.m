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
