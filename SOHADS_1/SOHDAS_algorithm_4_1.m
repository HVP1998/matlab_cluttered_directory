clc;clear;close all;
% ******************注释*************************
% 本程序旨在验证求解带有参数的方程或方程组
% ******************注释*************************
% syms a b c x y;
% sum=a*x^2+b*x+c;
% x=solve(sum==0)
syms x k U beta c1 c2 c3 c4 c5 c6 c7 c8 c9 c10;
A1con_3_a=U*exp(1i*k*x);
A1con_3_b=-(k^3*U^3*beta^2/16)*(k*x^2/2-1i*x)*exp(1i*k*x);  
A2_3=(k^2*U^2*beta*x/8)*exp(-2i*k*x);    
A2con_3=(k^2*U^2*beta*x/8)*exp(2i*k*x);
A3_3=(k^4*U^3*beta^2*x^2/32+1i*k^3*U^3*beta^2*x/24)*exp(-3i*k*x);
A1con_3_a_x=diff(A1con_3_a,x,1);
A1con_3_a_xx=diff(A1con_3_a,x,2);
A1con_3_b_x=diff(A1con_3_b,x,1);
A1con_3_b_xx=diff(A1con_3_b,x,2);
A2_3_x=diff(A2_3,x,1);
A2_3_xx=diff(A2_3,x,2);
A2con_3_x=diff(A2con_3,x,1);
A2con_3_xx=diff(A2con_3,x,2);
A3_3_x=diff(A3_3,x,1);
A3_3_xx=diff(A3_3,x,2);
%% 
% Eq1
sum1_right=(beta/2)*simplify(A2_3_xx*A1con_3_a_x+A2_3_x*A1con_3_a_xx);
param1_right=coeffs(sum1_right/(exp(-1i*k*x)),x);
A1_4_1=(c1*x^2+c2*x)*exp(-1i*k*x);
A1_4_1_xx=diff(A1_4_1,x,2);
sum1_left=-simplify(k^2*A1_4_1+A1_4_1_xx);
param1_left=coeffs(expand(sum1_left/exp(-1i*k*x)),x);
[c1,c2]=solve([param1_left(1)==param1_right(1),param1_left(2)==param1_right(2)],[c1,c2]);
A1_4_1=simplify(subs(A1_4_1));
%% 
% Eq2
sum2_right=(beta/2)*simplify(A1con_3_b_xx*A2_3_x+A1con_3_b_x*A2_3_xx+A2con_3_xx*A3_3_x+A2con_3_x*A3_3_xx);
param2_right=coeffs(sum2_right/(exp(-1i*k*x)),x);
A1_4_2=(c1*x^4+c2*x^3+c3*x^2+c4*x^1)*exp(-1i*k*x);
A1_4_2_xx=diff(A1_4_2,x,2);
sum2_left=-(k^2*A1_4_2+A1_4_2_xx);
param2_left=coeffs(expand(sum2_left/exp(-1i*k*x)),x);
[c1,c2,c3,c4]=solve([param2_left(1)==param2_right(1),param2_left(2)==param2_right(2),param2_left(3)==param2_right(3),param2_left(4)==param2_right(4)],[c1,c2,c3,c4]);
A1_4_2=simplify(subs(A1_4_2));