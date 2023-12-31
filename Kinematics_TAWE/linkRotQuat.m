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
