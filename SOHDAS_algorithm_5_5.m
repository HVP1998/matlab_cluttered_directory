clc;clear;close all;
% ******************注释*************************
% 本程序旨在验证求解带有参数的方程或方程组
% ******************注释*************************
% syms a b c x y;
% sum=a*x^2+b*x+c;
% x=solve(sum==0)
syms x k U beta;
syms c1_1 c1_2 c1_3 c1_4;
syms c2_1 c2_2 c2_3 c2_4 c2_5 c2_6;
syms c3_1 c3_2 c3_3 c3_4 c3_5 c3_6 c3_7 c3_8;
syms c4_1 c4_2 c4_3 c4_4 c4_5 c4_6 c4_7 c4_8 c4_9 c4_10;
A1_4_a=U*exp(-1i*k*x);
A1_4_b=-(k^3*U^3*beta^2/16)*(k*x^2/2+1i*x)*exp(-1i*k*x);
A1_4_c=(k^5*U^5*beta^4/128)*(-k^3*x^4/8-11*1i*k^2*x^3/12+5*k*x^2/12-5i*x/4)*exp(-1i*k*x);
A1con_4_a=U*exp(1i*k*x);
A1con_4_b=-(k^3*U^3*beta^2/16)*(k*x^2/2-1i*x)*exp(1i*k*x);
A1con_4_c=(k^5*U^5*beta^4/128)*(-k^3*x^4/8+11*1i*k^2*x^3/12+5*k*x^2/12+5i*x/4)*exp(1i*k*x);
A2_4_a=k^2*U^2*beta*x/8*exp(-2i*k*x);                         
A2_4_b=-(k^4*U^4*beta^3/24)*(k^2*x^3/4+9*1i*k*x^2/8-x)*exp(-2i*k*x);                      
A2_4_c=(k^6*U^6*beta^5/8192)*(7*k^4*x^5/5+3i*k^3*x^4/4+259*k^2*x^3/36-643i*k*x^2/48+1661*x/96)*exp(-2i*k*x);
A2con_4_a=(k^2*U^2*beta*x/8)*exp(2i*k*x);
A2con_4_b=(k^4*U^4*beta^3/24)*(-k^2*x^3/4+9i*k*x^2/8+x)*exp(2i*k*x);
A2con_4_c=(k^6*U^6*beta^5/8192)*(7*k^4*x^5/5-3i*k^3*x^4/4+259*k^2*x^3/36+643*1i*k*x^2/48+1661*x/96)*exp(2i*k*x);
A3_4_a=(k^3*U^3*beta^2/8)*(k*x^2/4+1i*x/3)*exp(-3i*k*x);                    
A3_4_b=(k^5*U^5*beta^4/64)*(-k^3*x^4/32-5*1i*k^2*x^3/25+k*x^2/3+13*1i*x/72)*exp(-3i*k*x);
A4_4_a=(k^4*U^4*beta^3/32)*(k^2*x^3/2+1i*k*x^2-7*x/12)*exp(-4i*k*x);                    
A4_4_b=(k^6*U^6*beta^5/4096)*(-3*k^4*x^5/5-39i*k^3*x^4/8+1801*k^2*x^3/144+4535i*k*x^2/384-6217*x/1536)*exp(-4i*k*x);
A1_4_a_x=diff(A1_4_a,x,1);
A1_4_a_xx=diff(A1_4_a,x,2);
A1_4_b_x=diff(A1_4_b,x,1);
A1_4_b_xx=diff(A1_4_b,x,2);
A1_4_c_x=diff(A1_4_c,x,1);
A1_4_c_xx=diff(A1_4_c,x,2);
A1con_4_a_x=diff(A1con_4_a,x,1);
A1con_4_a_xx=diff(A1con_4_a,x,2);
A1con_4_b_x=diff(A1con_4_b,x,1);
A1con_4_b_xx=diff(A1con_4_b,x,2);
A1con_4_c_x=diff(A1con_4_c,x,1);
A1con_4_c_xx=diff(A1con_4_c,x,2);
A2_4_a_x=diff(A2_4_a,x,1);
A2_4_a_xx=diff(A2_4_a,x,2);
A2_4_b_x=diff(A2_4_b,x,1);
A2_4_b_xx=diff(A2_4_b,x,2);
A2_4_c_x=diff(A2_4_c,x,1);
A2_4_c_xx=diff(A2_4_c,x,2);
A2con_4_a_x=diff(A2con_4_a,x,1);
A2con_4_a_xx=diff(A2con_4_a,x,2);
A2con_4_b_x=diff(A2con_4_b,x,1);
A2con_4_b_xx=diff(A2con_4_b,x,2);
A2con_4_c_x=diff(A2con_4_c,x,1);
A2con_4_c_xx=diff(A2con_4_c,x,2);
A3_4_a_x=diff(A3_4_a,x,1);
A3_4_a_xx=diff(A3_4_a,x,2);
A3_4_b_x=diff(A3_4_b,x,1);
A3_4_b_xx=diff(A3_4_b,x,2);
A4_4_a_x=diff(A4_4_a,x,1);
A4_4_a_xx=diff(A4_4_a,x,2);
A4_4_b_x=diff(A4_4_b,x,1);
A4_4_b_xx=diff(A4_4_b,x,2);
%%
% Eq1
sum1_right=(beta/2)*simplify(A1_4_a_x*A4_4_a_xx+A1_4_a_xx*A4_4_a_x+A2_4_a_x*A3_4_a_xx+A2_4_a_xx*A3_4_a_x);
param1_right=coeffs(sum1_right/(exp(-5i*k*x)),x);
A5_5_1=(c1_1*x^4+c1_2*x^3+c1_3*x^2+c1_4*x^1)*exp(-5i*k*x);
A5_5_1_xx=diff(A5_5_1,x,2);
sum1_left=-simplify(25*k^2*A5_5_1+A5_5_1_xx);
param1_left=coeffs(expand(sum1_left/exp(-5i*k*x)),x);
[c1_1,c1_2,c1_3,c1_4]=solve([param1_left(1)==param1_right(1),param1_left(2)==param1_right(2),param1_left(3)==param1_right(3),param1_left(4)==param1_right(4)],[c1_1,c1_2,c1_3,c1_4]);
A5_5_1=simplify(subs(A5_5_1));
%% 
% Eq2
sum2_right=(beta/2)*simplify(A1_4_a_xx*A4_4_b_x+A1_4_a_x*A4_4_b_xx+A1_4_b_xx*A4_4_a_x+A1_4_b_x*A4_4_a_xx+A2_4_a_xx*A3_4_b_x+A2_4_a_x*A3_4_b_xx+A2_4_b_xx*A3_4_a_x+A2_4_b_x*A3_4_a_xx);
param2_right=coeffs(sum2_right/(exp(-5i*k*x)),x);
A5_5_2=(c2_1*x^6+c2_2*x^5+c2_3*x^4+c2_4*x^3+c2_5*x^2+c2_6*x^1)*exp(-5i*k*x);
A5_5_2_xx=diff(A5_5_2,x,2);
sum2_left=-simplify(25*k^2*A5_5_2+A5_5_2_xx);
param2_left=coeffs(expand(sum2_left/exp(-5i*k*x)),x);
C2=solve([param2_left(1)==param2_right(1),param2_left(2)==param2_right(2),param2_left(3)==param2_right(3),param2_left(4)==param2_right(4),param2_left(5)==param2_right(5),param2_left(6)==param2_right(6)],[c2_1,c2_2,c2_3,c2_4,c2_5,c2_6]);
A5_5_2=simplify(subs(A5_5_2,C2));
%% 
% Eq3
sum3_right=(beta/2)*simplify(A1_4_b_xx*A4_4_b_x+A1_4_b_x*A4_4_b_xx+A1_4_c_xx*A4_4_a_x+A1_4_c_x*A4_4_a_xx+A2_4_b_xx*A3_4_b_x+A2_4_b_x*A3_4_b_xx+A2_4_c_xx*A3_4_a_x+A2_4_c_x*A3_4_a_xx);
param3_right=coeffs(sum3_right/(exp(-5i*k*x)),x);
A5_5_3=(c3_1*x^8+c3_2*x^7+c3_3*x^6+c3_4*x^5+c3_5*x^4+c3_6*x^3+c3_7*x^2+c3_8*x^1)*exp(-5i*k*x);
A5_5_3_xx=diff(A5_5_3,x,2);
sum3_left=-simplify(25*k^2*A5_5_3+A5_5_3_xx);
param3_left=coeffs(expand(sum3_left/exp(-5i*k*x)),x);
[c3_1,c3_2,c3_3,c3_4,c3_5,c3_6,c3_7,c3_8]=solve([param3_left(1)==param3_right(1),param3_left(2)==param3_right(2),param3_left(3)==param3_right(3),param3_left(4)==param3_right(4),param3_left(5)==param3_right(5),param3_left(6)==param3_right(6),param3_left(7)==param3_right(7),param3_left(8)==param3_right(8)],[c3_1,c3_2,c3_3,c3_4,c3_5,c3_6,c3_7,c3_8]);
A5_5_3=simplify(subs(A5_5_3));
%% 
% Eq4
sum4_right=(beta/2)*simplify(A1_4_c_x*A4_4_b_xx+A1_4_c_xx*A4_4_b_x+A2_4_c_x*A3_4_b_xx+A2_4_c_xx*A3_4_b_x);
param4_right=coeffs(sum4_right/(exp(-5i*k*x)),x);
A5_5_4=(c4_1*x^10+c4_2*x^9+c4_3*x^8+c4_4*x^7+c4_5*x^6+c4_6*x^5+c4_7*x^4+c4_8*x^3+c4_9*x^2+c4_10*x^1)*exp(-5i*k*x);
A5_5_4_xx=diff( A5_5_4,x,2);
sum4_left=-simplify(25*k^2*A5_5_4+A5_5_4_xx);
param4_left=coeffs(expand(sum4_left/exp(-5i*k*x)),x);
[c4_1,c4_2,c4_3,c4_4,c4_5,c4_6,c4_7,c4_8,c4_9,c4_10]=solve([param4_left(1)==param4_right(1),param4_left(2)==param4_right(2),param4_left(3)==param4_right(3),param4_left(4)==param4_right(4),param4_left(5)==param4_right(5),param4_left(6)==param4_right(6),param4_left(7)==param4_right(7),param4_left(8)==param4_right(8),param4_left(9)==param4_right(9),param4_left(10)==param4_right(10)],[c4_1,c4_2,c4_3,c4_4,c4_5,c4_6,c4_7,c4_8,c4_9,c4_10]);
A5_5_4=simplify(subs(A5_5_4));