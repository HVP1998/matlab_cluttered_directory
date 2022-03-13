% *****************************************************
% �����µ���ⷽ����߽ײ����㶯��
% *****************************************************
clc;clear;close all;
f=2e6;                                                                 %���䳬��Ƶ��f=1MHz
omega=2*pi*f;                                                    %���ٶȦ�=2*��*f
c=1491;                                                               %ˮ�е�����c=1491m/s
k=omega/c;                                                         %����k=��/c
t=2;                                                                     %ʱ��ȡ�����������2s
x=0:0.001:0.15;                                                    %��������ȡ0��0.15m
beta=6.88;                                                            %ˮ�еķ����Բ�����=3.5
U=87.5e-9;                                                           %��Դ���𶯷���U=100nm
theta=omega*t-k*x;
y2_2=(k^2*U^2*beta/8)*x.*cos(2*theta);
y2_2_amplitude_cos=(k^2*U^2*beta/8)*x;
y2_2_amplitude_sin=0;
y2_2_amplitude=((y2_2_amplitude_cos).^2+(y2_2_amplitude_sin).^2).^0.5;
y2_3=y2_2;
y2_3_amplitude_cos=y2_2_amplitude_cos;
y2_3_amplitude_sin=y2_2_amplitude_sin;
y2_3_amplitude=((y2_3_amplitude_cos).^2+(y2_3_amplitude_sin).^2).^0.5;
y2_4_a=y2_2;
y2_4_a_amplitude_cos=(k^2*U^2*beta/8)*x;
y2_4_a_amplitude_sin=0;
y2_4_b=(k^4*U^4*beta^3/24)*(-(k^2/4)*(x.^3).*cos(2*theta)+(9*k/8)*(x.^2).*sin(2*theta)+x.*cos(2*theta));
y2_4_b_amplitude_cos=(k^4*U^4*beta^3/24)*(-(k^2/4)*(x.^3)+x);
y2_4_b_amplitude_sin=(k^4*U^4*beta^3/24)*(9*k/8)*(x.^2);
y2_4_c=(k^6*U^6*beta^5/8192)*((7*k^4/5)*(x.^5).*cos(2*theta)-(3*k^3/4)*(x.^4).*sin(2*theta)+(259*k^2/36)*(x.^3).*cos(2*theta)+(643*k/48)*(x.^2).*sin(2*theta)+(1661/96)*x.*cos(2*theta));
y2_4_c_amplitude_cos=(k^6*U^6*beta^5/8192)*((7*k^4/5)*(x.^5)+(259*k^2/36)*(x.^3)+(1661/96)*x);
y2_4_c_amplitude_sin=(k^6*U^6*beta^5/8192)*(-(3*k^3/4)*(x.^4)+(643*k/48)*(x.^2));
y2_4=y2_4_a+y2_4_b+y2_4_c;
y2_4_amplitude_cos=y2_4_a_amplitude_cos+y2_4_b_amplitude_cos+y2_4_c_amplitude_cos;
y2_4_amplitude_sin=y2_4_a_amplitude_sin+y2_4_b_amplitude_sin+y2_4_c_amplitude_sin;
y2_4_amplitude=((y2_4_amplitude_cos).^2+(y2_4_amplitude_sin).^2).^0.5;
y2_5_a=y2_2;
y2_5_a_amplitude_cos=y2_2_amplitude_cos;
y2_5_a_amplitude_sin=y2_2_amplitude_sin;
y2_5_b=y2_4_b;
y2_5_b_amplitude_cos=y2_4_b_amplitude_cos;
y2_5_b_amplitude_sin=y2_4_b_amplitude_sin;
y2_5_c=-(U^6*beta^5*k^6/2949120)*(1728*k^4*(x.^5).*cos(2*theta)-8952*k^3*(x.^4).*sin(2*theta) - 9208*k^2*(x.^3).*cos(2*theta)-18042*k*(x.^2).*sin(2*theta) - 23299*x.*cos(2*theta));
y2_5_c_amplitude_cos=-(U^6*beta^5*k^6/2949120)*(1728*k^4*(x.^5) - 9208*k^2*(x.^3)- 23299*x);
y2_5_c_amplitude_sin=-(U^6*beta^5*k^6/2949120)*(-8952*k^3*(x.^4)-18042*k*(x.^2));
y2_5_d=(U^8*beta^7*k^8/211392921600)*(11512320*k^6*(x.^7).*cos(2*theta)+k^5*(x.^6)*10586240.*sin(2*theta)+ 86812992*k^4*(x.^5).*cos(2*theta)+k^3*(x.^4)*403709600.*sin(2*theta)+ 532493360*k^2*(x.^3).*cos(2*theta)+k*(x.^2)*153014750.*sin(2*theta) + 226904125*x.*cos(2*theta));
y2_5_d_amplitude_cos=(U^8*beta^7*k^8/211392921600)*(11512320*k^6*(x.^7)+ 86812992*k^4*(x.^5)+ 532493360*k^2*(x.^3)+ 226904125*x);
y2_5_d_amplitude_sin=(U^8*beta^7*k^8/211392921600)*(k^5*(x.^6)*10586240+k^3*(x.^4)*403709600+k*(x.^2)*153014750);
y2_5_e=-(U^10*beta^9*k^10/30440580710400)*(29191680*k^8*(x.^9).*cos(2*theta)-k^7*(x.^8)*83885760.*sin(2*theta)+745426560*k^6*(x.^7).*cos(2*theta)-k^5*(x.^6)*2746843680.*sin(2*theta)-109166736*k^4*(x.^5).*cos(2*theta)-k^3*(x.^4)*9346975680.*sin(2*theta)- 7403656120*k^2*(x.^3).*cos(2*theta)-k*(x.^2)*2079789390.*sin(2*theta)- 2370095805*x.*cos(2*theta));
y2_5_e_amplitude_cos=-(U^10*beta^9*k^10/30440580710400)*(29191680*k^8*(x.^9)+745426560*k^6*(x.^7)-109166736*k^4*(x.^5)- 7403656120*k^2*(x.^3)- 2370095805*x);
y2_5_e_amplitude_sin=-(U^10*beta^9*k^10/30440580710400)*(-k^7*(x.^8)*83885760-k^5*(x.^6)*2746843680-k^3*(x.^4)*9346975680-k*(x.^2)*2079789390);
y2_5_f=-(U^12*beta^11*k^12/685765402243891200)*(-3121348608*k^10*(x.^11).*cos(2*theta)+k^9*(x.^10)*30349541376.*sin(2*theta)+19383992320*k^8*(x.^9).*cos(2*theta)+k^7*(x.^8)*502196808960.*sin(2*theta)+1586121320960*k^6*(x.^7).*cos(2*theta)-k^5*(x.^6)*993841242240.*sin(2*theta)+ 3701692764000*k^4*(x.^5).*cos(2*theta)-k^3*(x.^4)*7696444371000.*sin(2*theta)- 3506369328000*k^2*(x.^3).*cos(2*theta)-k*(x.^2)*1310910158250.*sin(2*theta)- 827161317675*x.*cos(2*theta));
y2_5_f_amplitude_cos=-(U^12*beta^11*k^12/685765402243891200)*(-3121348608*k^10*(x.^11)+19383992320*k^8*(x.^9)+1586121320960*k^6*(x.^7)+ 3701692764000*k^4*(x.^5)- 3506369328000*k^2*(x.^3)- 827161317675*x);
y2_5_f_amplitude_sin=-(U^12*beta^11*k^12/685765402243891200)*(k^9*(x.^10)*30349541376+k^7*(x.^8)*502196808960-k^5*(x.^6)*993841242240-k^3*(x.^4)*7696444371000-k*(x.^2)*1310910158250);
y2_5=y2_5_a+y2_5_b+y2_5_c+y2_5_d+y2_5_e+y2_5_f;
y2_5_amplitude_cos=y2_5_a_amplitude_cos+y2_5_b_amplitude_cos+y2_5_c_amplitude_cos+y2_5_d_amplitude_cos+y2_5_e_amplitude_cos+y2_5_f_amplitude_cos;
y2_5_amplitude_sin=y2_5_a_amplitude_sin+y2_5_b_amplitude_sin+y2_5_c_amplitude_sin+y2_5_d_amplitude_sin+y2_5_e_amplitude_sin+y2_5_f_amplitude_sin;
y2_5_amplitude=((y2_5_amplitude_cos).^2+(y2_5_amplitude_sin).^2).^0.5;
% figure(1);
% plot(x,y2_2);
% hold on;
% plot(x,y2_2);
% legend('old','new')
% title('�ڶ��ε�������λ��')
% figure(2);
% plot(x,y2_3);
% hold on;
% plot(x,y2_3);
% legend('old','new')
% title('�����ε�������λ��')
% figure(3);
% plot(x,y2_4_a+y2_4_b);
% hold on;
% plot(x,y2_4);
% legend('old','new')
% title('���Ĵε�������λ��')
% figure(4);
% plot(x,y2_4_a+y2_4_b);
% hold on;
% plot(x,y2_5);
% legend('old','new')
% title('����ε�������λ��')
figure(5);
plot(x,y2_2_amplitude/U);
hold on;
plot(x,y2_3_amplitude/U);
hold on;
plot(x,y2_4_amplitude/U);
hold on;
plot(x,y2_5_amplitude/U);
legend("���ε���","���ε���","�Ĵε���","��ε���");
title('����λ�Ʒ�ֵ')