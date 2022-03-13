% *****************************************************
% 复刻文献
% *****************************************************
clc;clear;close all;
f=2e6;                                                                 %发射超声频率f=2MHz
omega=2*pi*f;                                                    %角速度ω=2*π*f
c=1491;                                                               %水中的声速c=1491m/s
k=omega/c;                                                         %波数k=ω/c
t=2;                                                                     %时间取声波传播后的2s
x=0:0.001:0.15;                                                    %传播距离取0到0.15m
beta=6.88;                                                           %水中的非线性参数β=6.88
U=100e-9;                                                           %声源的震动幅度U=100nm
theta=omega*t-k*x;
% % 第一次迭代
% y1_1=U*cos(theta);
% % 第二次迭代
% y1_2=y1_1;
% y2_2=(k^2*U^2*beta/8)*x.*cos(2*theta);
% %第三次迭代
% y1_3_a=(-k^3*U^3*beta^2/16)*((k/2)*(x.^2).*cos(theta)-x.*sin(theta));
% y1_3_b=y1_1;
% y1_3=y1_3_a+y1_3_b;
% y2_3=y2_2;
% y3_3=(k^3*U^3*beta^2/16)*((k/2)*(x.^2).*cos(3*theta)-(2/3)*x.*sin(3*theta));
% %第四次迭代
% y1_4_b=y1_3_a;
% y1_4_c=(k^5*U^5*beta^4/128)*(-(k^3/8)*(x.^4).*cos(theta)+(11*k^2/12)*(x.^3).*sin(theta)-(25*k/18)*(x.^2).*cos(theta)-(5/9)*x.*sin(theta));
% y1_4_a=y1_1;
% y1_4=y1_4_a+y1_4_b+y1_4_c;
% y2_4_a=y2_2;
% y2_4_b=(k^4*U^4*beta^3/24)*(-(k^2/4)*(x.^3).*cos(2*theta)+(9*k/8)*(x.^2).*sin(2*theta)+x.*cos(2*theta));
% y2_4_c=(k^6*U^6*beta^5/8192)*((7*k^4/5)*(x.^5).*cos(2*theta)+(17*k^3/4)*(x.^4).*sin(2*theta)+(79*k^2/36)*(x.^3).*cos(2*theta)+(463*k/48)*(x.^2).*sin(2*theta)+(1841/96)*x.*cos(2*theta));
% y2_4=y2_4_a+y2_4_b+y2_4_c;
% y3_4_a=y3_3;
% y3_4_b=(k^5*U^5*beta^4/64)*((-k^3/32)*(x.^4).*cos(3*theta)+(5*k^2/24)*(x.^3).*sin(3*theta)+(k/3)*(x.^2).*cos(3*theta)-(13/72)*x.*sin(3*theta));
% y3_4=y3_4_a+y3_4_b;
% y4_4_a=(k^4*U^4*beta^3/96)*((k^2)*(x.^3).*cos(4*theta)-(3*k)*(x.^2).*sin(4*theta)-(7/4)*x.*cos(4*theta));
% y4_4_b=(k^6*U^6*beta^5/4096)*(-(3*k^4/5)*(x.^5).*cos(4*theta)+(39*k^3/8)*(x.^4).*sin(4*theta)+(1801*k^2/144)*(x.^3).*cos(4*theta)-(4535*k/384)*(x.^2).*sin(4*theta)-(6217/192)*x.*cos(4*theta));
% y4_4=y4_4_a+y4_4_b;
u2_2=(k^2*U^2*beta/8)*x.*cos(2*theta);
u2_2_amplitude_cos=(k^2*U^2*beta/8)*x;
u2_2_amplitude_sin=0;
u2_2_amplitude=((u2_2_amplitude_cos).^2+(u2_2_amplitude_sin).^2).^0.5;
u2_4_a=u2_2;
u2_4_b=-(k^6*U^4*beta^3/96)*(x.^3).*cos(2*theta)+(3*k^5*U^4*beta^3/64)*(x.^2).*sin(2*theta)+(k^4*U^4*beta^3/96)*x.*cos(2*theta);
u2_4=u2_4_a+u2_4_b;
u2_4_amplitude_cos=u2_2_amplitude_cos-(k^6*U^4*beta^3/96)*(x.^3)+(k^4*U^4*beta^3/96)*x;
u2_4_amplitude_sin=u2_2_amplitude_sin+(3*k^5*U^4*beta^3/64)*(x.^2);
u2_4_amplitude=((u2_4_amplitude_cos).^2+(u2_4_amplitude_sin).^2).^0.5;
u2_6_a=u2_2;
u2_6_b=u2_4_b;
u2_6_c=(k^10*U^6*beta^5*(x.^5)/3072-85*k^8*U^6*beta^5*(x.^3)/3072+175*k^6*U^6*beta^5*x/6144).*cos(2*theta)+(-k^9*U^6*beta^5*(x.^4)/192+325*k^7*U^6*beta^5*(x.^2)/6144).*sin(2*theta);
u2_6=u2_6_a+u2_6_b+u2_6_c;
u2_6_amplitude_cos=u2_4_amplitude_cos+(k^10*U^6*beta^5*(x.^5)/3072-85*k^8*U^6*beta^5*(x.^3)/3072+175*k^6*U^6*beta^5*x/6144);
u2_6_amplitude_sin=u2_4_amplitude_sin+(-k^9*U^6*beta^5*(x.^4)/192+325*k^7*U^6*beta^5*(x.^2)/6144);
u2_6_amplitude=((u2_6_amplitude_cos).^2+(u2_6_amplitude_sin).^2).^0.5;
u2_8_a=u2_2;
u2_8_b=u2_4_b;
u2_8_c=u2_6_c;
u2_8_d=(-k^14*U^8*beta^7*(x.^7)/184320+91*k^12*U^8*beta^7*(x.^5)/36864+6935*k^10*U^8*beta^7*(x.^3)/147456++79*k^8*U^8*beta^7*x/3072).*cos(2*theta)+(23*k^13*U^8*beta^7*(x.^6)/122880+569*k^11*U^8*beta^7*(x.^4)/36864+18397*k^9*U^8*beta^7*(x.^2)/294912).*sin(2*theta);
u2_8=u2_8_a+u2_8_b+u2_8_c+u2_8_d;
u2_8_amplitude_cos=u2_6_amplitude_cos+(-k^14*U^8*beta^7*(x.^7)/184320+91*k^12*U^8*beta^7*(x.^5)/36864+6935*k^10*U^8*beta^7*(x.^3)/147456++79*k^8*U^8*beta^7*x/3072);
u2_8_amplitude_sin=u2_6_amplitude_sin+(23*k^13*U^8*beta^7*(x.^6)/122880+569*k^11*U^8*beta^7*(x.^4)/36864+18397*k^9*U^8*beta^7*(x.^2)/294912);
u2_8_amplitude=((u2_8_amplitude_cos).^2+(u2_8_amplitude_sin).^2).^0.5;
% figure(1);
% plot(x,u2_2);
% title('第二次迭代二阶位移')
% figure(2);
% plot(x,u2_4);
% title('第四次迭代二阶位移')
% figure(3);
% plot(x,u2_6);
% title('第六次迭代二阶位移')
% figure(4);
% plot(x,u2_8);
% title('第八次迭代二阶位移')
figure(5);
plot(x,u2_2_amplitude/U);
hold on;
plot(x,u2_4_amplitude/U);
hold on;
plot(x,u2_6_amplitude/U);
hold on;
plot(x,u2_8_amplitude/U);
legend("二次迭代","四次迭代","六次迭代","八次迭代");
title('二阶位移幅值')