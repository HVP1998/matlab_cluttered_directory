clc;clear
f=2e6;                                                                 %���䳬��Ƶ��f=1MHz
omega=2*pi*f;                                                    %���ٶȦ�=2*��*f
c=1491;                                                               %ˮ�е�����c=1491m/s
k=omega/c;                                                         %����k=��/c
t=2;                                                                     %ʱ��ȡ�����������2s
x=0:0.001:0.15;                                                    %��������ȡ0��0.15m
beta=3.5;                                                            %ˮ�еķ����Բ�����=3.5
U=100e-9;                                                           %��Դ���𶯷���U=100nm
theta=omega*t-k*x;
y1_1=U*cos(theta);
y1_2=y1_1;y2_2=0;
y1_3_a=(-k^3*U^3*beta^2/16)*((k/2)*(x.^2).*cos(theta)-x.*sin(theta));
y1_3_b=(k^5*U^5*beta^4/128)*(-(k^3/8)*(x.^4).*cos(theta)+(11*k^2/12)*(x.^3).*sin(theta)-(25*k/18)*(x.^2).*cos(theta)-(5/9)*x.*sin(theta));
y2_4_a=(k^4*U^4*beta^3/24)*(-(k^2/4)*(x.^3).*cos(2*theta)+(9*k/8)*(x.^2).*sin(2*theta)+x.*cos(2*theta));