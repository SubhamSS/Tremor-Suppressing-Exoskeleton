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