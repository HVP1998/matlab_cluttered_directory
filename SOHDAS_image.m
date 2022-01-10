clc;clear
f=2e6;                                                                 %发射超声频率f=1MHz
omega=2*pi*f;                                                    %角速度ω=2*π*f
c=1491;                                                               %水中的声速c=1491m/s
k=omega/c;                                                         %波数k=ω/c
t=2;                                                                     %时间取声波传播后的2s
x=0:0.001:0.15;                                                    %传播距离取0到0.15m
beta=3.5;                                                            %水中的非线性参数β=3.5
U=100e-9;                                                           %声源的震动幅度U=100nm
theta=omega*t-k*x;
y1_1=U*cos(theta);
y1_2=y1_1;y2_2=0;
y1_3_a=(-k^3*U^3*beta^2/16)*((k/2)*(x.^2).*cos(theta)-x.*sin(theta));
y1_3_b=(k^5*U^5*beta^4/128)*(-(k^3/8)*(x.^4).*cos(theta)+(11*k^2/12)*(x.^3).*sin(theta)-(25*k/18)*(x.^2).*cos(theta)-(5/9)*x.*sin(theta));
y2_4_a=(k^4*U^4*beta^3/24)*(-(k^2/4)*(x.^3).*cos(2*theta)+(9*k/8)*(x.^2).*sin(2*theta)+x.*cos(2*theta));