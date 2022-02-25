clc;clear;close all;
% ******************注释*************************
% 本程序旨在验证求解带有参数的方程或方程组
% ******************注释*************************
% syms a b c x y;
% sum=a*x^2+b*x+c;
% x=solve(sum==0)
syms x k U beta;
syms c1_1 c1_2 c1_3 c1_4 c1_5 c1_6;
syms c2_1 c2_2 c2_3 c2_4 c2_5 c2_6 c2_7 c2_8;
syms c3_1 c3_2 c3_3 c3_4 c3_5 c3_6 c3_7 c3_8 c3_9 c3_10;
A3con_4_a=(k^3*U^3*beta^2/8)*(k*x^2/4-1i*x/3)*exp(3i*k*x);                    
A3con_4_b=(k^5*U^5*beta^4/64)*(-k^3*x^4/32+5*1i*k^2*x^3/25+k*x^2/3-13*1i*x/72)*exp(3i*k*x);
A4_4_a=(k^4*U^4*beta^3/32)*(k^2*x^3/2+1i*k*x^2-7*x/12)*exp(-4i*k*x);                    
A4_4_b=(k^6*U^6*beta^5/4096)*(-3*k^4*x^5/5-39i*k^3*x^4/8+1801*k^2*x^3/144+4535i*k*x^2/384-6217*x/1536)*exp(-4i*k*x);
A3con_4_a_x=diff(A3con_4_a,x,1);
A3con_4_a_xx=diff(A3con_4_a,x,2);
A3con_4_b_x=diff(A3con_4_b,x,1);
A3con_4_b_xx=diff(A3con_4_b,x,2);
A4_4_a_x=diff(A4_4_a,x,1);
A4_4_a_xx=diff(A4_4_a,x,2);
A4_4_b_x=diff(A4_4_b,x,1);
A4_4_b_xx=diff(A4_4_b,x,2);
%% 
% Eq1
sum1_right=(beta/2)*simplify(A3con_4_a_x*A4_4_a_xx+A3con_4_a_xx*A4_4_a_x);
param1_right=coeffs(sum1_right/exp(-k*x*1i),x);
A1_5_1=(c1_1*x^6+c1_2*x^5+c1_3*x^4+c1_4*x^3+c1_5*x^2+c1_6*x)*exp(-1i*k*x);
A1_5_1_xx=diff(A1_5_1,x,2);
sum1_left=-(k^2*A1_5_1+A1_5_1_xx);
param1_left=coeffs(simplify(sum1_left/exp(-k*x*1i)),x);
C1=solve([param1_left(1)==param1_right(1),param1_left(2)==param1_right(2),param1_left(3)==param1_right(3),param1_left(4)==param1_right(4),param1_left(5)==param1_right(5),param1_left(6)==param1_right(6)],[c1_1,c1_2,c1_3,c1_4,c1_5,c1_6]);
A1_5_1=simplify(subs(A1_5_1,C1));
%% 
% Eq2 
sum2_right=(beta/2)*simplify(A3con_4_a_x*A4_4_b_xx+A3con_4_a_xx*A4_4_b_x+A3con_4_b_x*A4_4_a_xx+A3con_4_b_xx*A4_4_a_x);
param2_right=coeffs(sum2_right/exp(-k*x*1i),x);
A1_5_2=(c2_1*x^8+c2_2*x^7+c2_3*x^6+c2_4*x^5+c2_5*x^4+c2_6*x^3+c2_7*x^2+c2_8*x)*exp(-1i*k*x);
A1_5_2_xx=diff(A1_5_2,x,2);
sum2_left=-(k^2*A1_5_2+A1_5_2_xx);
param2_left=coeffs(simplify(sum2_left/exp(-k*x*1i)),x);
C2=solve([param2_left(1)==param2_right(1),param2_left(2)==param2_right(2),param2_left(3)==param2_right(3),param2_left(4)==param2_right(4),param2_left(5)==param2_right(5),param2_left(6)==param2_right(6),param2_left(7)==param2_right(7),param2_left(8)==param2_right(8)],[c2_1,c2_2,c2_3,c2_4,c2_5,c2_6,c2_7,c2_8]);
A1_5_2=simplify(subs(A1_5_2,C2));
%% 
%Eq3
sum3_right=(beta/2)*simplify(A3con_4_b_x*A4_4_b_xx+A3con_4_b_xx*A4_4_b_x);
param3_right=coeffs(sum3_right/exp(-k*x*1i),x);
A1_5_3=(c3_1*x^10+c3_2*x^9+c3_3*x^8+c3_4*x^7+c3_5*x^6+c3_6*x^5+c3_7*x^4+c3_8*x^3+c3_9*x^2+c3_10*x^1)*exp(-1i*k*x);
A1_5_3_xx=diff(A1_5_3,x,2);
sum3_left=-(k^2*A1_5_3+A1_5_3_xx);
param3_left=coeffs(simplify(sum3_left/exp(-k*x*1i)),x);
C3=solve([param3_left(1)==param3_right(1),param3_left(2)==param3_right(2),param3_left(3)==param3_right(3),param3_left(4)==param3_right(4),param3_left(5)==param3_right(5),param3_left(6)==param3_right(6),param3_left(7)==param3_right(7),param3_left(8)==param3_right(8),param3_left(9)==param3_right(9),param3_left(10)==param3_right(10)],[c3_1,c3_2,c3_3,c3_4,c3_5,c3_6,c3_7,c3_8,c3_9,c3_10]);
A1_5_3=simplify(subs(A1_5_3,C3));