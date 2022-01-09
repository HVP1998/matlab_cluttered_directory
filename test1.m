clc;clear;close all;
% ******************注释*************************
% 本程序旨在验证求解带有参数的方程或方程组
% ******************注释*************************
% syms a b c x y;
% sum=a*x^2+b*x+c;
% x=solve(sum==0)
syms x k U beta c1 c2 c3 c4;
A1con_4_a=U*exp(1i*k*x);
A1con_4_b=-(k^3*U^3*beta^2/16)*(k*x^2/2-1i*x)*exp(1i*k*x);
% coeffs(A1con_4_b/exp(1i*k*x),x)                       %求公式A1con_4_b中x的系数
A1con_4_c=(k^5*U^5*beta^4/128)*(-k^3*x^4/8+11*1i*k^2*x^3/12-25*k*x^2/18-5*1i/9)*exp(1i*k*x);
% coeffs(A1con_4_c/((k^5*U^5*beta^4/128)*exp(1i*k*x)),x)
A2_4_a=k^2*U^2*beta*x/8*exp(-2i*k*x);
% coeffs(A2_4_a/exp(-2i*k*x),x)                       
A2_4_b=-(k^4*U^4*beta^3/24)*(k^2*x^3/4+9*1i*k*x^2/8-x)*exp(-2i*k*x);
% coeffs(A2_4_b/((k^4*U^4*beta^3/24)*exp(-2i*k*x)),x)                      
A2_4_c=(k^6*U^6*beta^5/8192)*(9*k^4*x^5/5-17*1i*k^3*x^4/4+79*k^2*x^3/36-463*1i*k*x^2/48+1841*x/96)*exp(-2i*k*x);
% coeffs(A2_4_c/((k^6*U^6*beta^5/8192)*exp(-2i*k*x)),x)
A1con_4_a_x=diff(A1con_4_a,x,1);
A2_4_b_xx=diff(A2_4_b,x,2);
% coeffs(A2_4_b_xx*A1con_4_a_x/(exp(-1i*k*x)),x)
A1con_4_a_xx=diff(A1con_4_a,x,2);
A2_4_b_x=diff(A2_4_b,x,1);
% coeffs(A2_4_b_x*A1con_4_a_xx/(exp(-1i*k*x)),x)
A1con_4_b_x=diff(A1con_4_b,x,1);
A2_4_a_xx=diff(A2_4_a,x,2);
% coeffs(A2_4_a_xx*A1con_4_b_x/(exp(-1i*k*x)),x)
A1con_4_b_xx=diff(A1con_4_b,x,2);
A2_4_a_x=diff(A2_4_a,x,1);
% coeffs(A2_4_a_x*A1con_4_b_xx/(exp(-1i*k*x)),x)
sum2_right=(beta/2)*(A2_4_b_xx*A1con_4_a_x+A2_4_b_x*A1con_4_a_xx+A2_4_a_xx*A1con_4_b_x+A2_4_a_x*A1con_4_b_xx);
param2_right=coeffs(sum2_right/(exp(-1i*k*x)),x);
A1_5=(c1*x^4+c2*x^3+c3*x^2+c4*x)*exp(-1i*k*x);
A1_5_xx=diff(A1_5,x,2);
sum2_left=-(k^2*A1_5+A1_5_xx);
param2_left=coeffs(expand(sum2_left/exp(-1i*k*x)),x);
[c1,c2,c3,c4]=solve([param2_left(1)==param2_right(1),param2_left(2)==param2_right(2),param2_left(3)==param2_right(3),param2_left(4)==param2_right(4)],[c1,c2,c3,c4]);