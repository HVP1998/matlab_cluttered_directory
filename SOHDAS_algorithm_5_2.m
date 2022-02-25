clc;clear;close all;
% ******************注释*************************
% 本程序旨在验证求解带有参数的方程或方程组
% ******************注释*************************
% syms a b c x y;
% sum=a*x^2+b*x+c;
% x=solve(sum==0)
syms x k U beta c1 c2 c3 c4 c5 c6 c7 c8 c9 c10;
syms c1_1;
syms c2_1 c2_2 c2_3;
syms c3_1 c3_2 c3_3 c3_4 c3_5;
syms c4_1 c4_2 c4_3 c4_4 c4_5 c4_6 c4_7;
syms c5_1 c5_2 c5_3 c5_4 c5_5 c5_6 c5_7 c5_8 c5_9;
syms c6_1 c6_2 c6_3 c6_4 c6_5 c6_6 c6_7 c6_8 c6_9 c6_10 c6_11;
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
sum1_right=(beta/2)*simplify(A1_4_a_xx*A1_4_a_x);
param1_right=coeffs(sum1_right/(exp(-2i*k*x)),x);
A2_5_1=(c1_1*x)*exp(-2i*k*x);
A2_5_1_xx=diff(A2_5_1,x,2);
sum1_left=-simplify(4*k^2*A2_5_1+A2_5_1_xx);
param1_left=coeffs(expand(sum1_left/exp(-2i*k*x)),x);
c1_1=solve(param1_left(1)==param1_right(1),c1_1);
A2_5_1=simplify(subs(A2_5_1));
%% 
% Eq2
sum2_right=(beta/2)*simplify(A1_4_b_xx*A1_4_a_x+A1_4_b_x*A1_4_a_xx+A1con_4_a_xx*A3_4_a_x+A1con_4_a_x*A3_4_a_xx);
param2_right=coeffs(sum2_right/(exp(-2i*k*x)),x);
A2_5_2=(c2_1*x^3+c2_2*x^2+c2_3*x^1)*exp(-2i*k*x);
A2_5_2_xx=diff(A2_5_2,x,2);
sum2_left=-simplify(4*k^2*A2_5_2+A2_5_2_xx);
param2_left=coeffs(expand(sum2_left/exp(-2i*k*x)),x);
C2=solve([param2_left(1)==param2_right(1),param2_left(2)==param2_right(2),param2_left(3)==param2_right(3)],[c2_1,c2_2,c2_3]);
A2_5_2=simplify(subs(A2_5_2,C2));
%% 
% Eq3
sum3_right=(beta/2)*simplify(A1_4_b_x*A1_4_b_xx+A1_4_a_x*A1_4_c_xx+A1_4_a_xx*A1_4_c_x+A1con_4_a_x*A3_4_b_xx+A1con_4_a_xx*A3_4_b_x+A1con_4_b_x*A3_4_a_xx+A1con_4_b_xx*A3_4_a_x+A2con_4_a_x*A4_4_a_xx+A2con_4_a_xx*A4_4_a_x);
param3_right=coeffs(sum3_right/(exp(-2i*k*x)),x);
A2_5_3=(c3_1*x^5+c3_2*x^4+c3_3*x^3+c3_4*x^2+c3_5*x^1)*exp(-2i*k*x);
A2_5_3_xx=diff(A2_5_3,x,2);
sum3_left=-simplify(4*k^2*A2_5_3+A2_5_3_xx);
param3_left=coeffs(expand(sum3_left/exp(-2i*k*x)),x);
[c3_1,c3_2,c3_3,c3_4,c3_5]=solve([param3_left(1)==param3_right(1),param3_left(2)==param3_right(2),param3_left(3)==param3_right(3),param3_left(4)==param3_right(4),param3_left(5)==param3_right(5)],[c3_1,c3_2,c3_3,c3_4,c3_5]);
A2_5_3=simplify(subs(A2_5_3));
%% 
% Eq4
sum4_right=(beta/2)*simplify(A1_4_b_x*A1_4_c_xx+A1_4_b_xx*A1_4_c_x+A1con_4_b_x*A3_4_b_xx+A1con_4_b_xx*A3_4_b_x+A1con_4_c_x*A3_4_a_xx+A1con_4_c_xx*A3_4_a_x+A2con_4_a_x*A4_4_b_xx+A2con_4_a_xx*A4_4_b_x+A2con_4_b_x*A4_4_a_xx+A2con_4_b_xx*A4_4_a_x);
param4_right=coeffs(sum4_right/(exp(-2i*k*x)),x);
A2_5_4=(c4_1*x^7+c4_2*x^6+c4_3*x^5+c4_4*x^4+c4_5*x^3+c4_6*x^2+c4_7*x^1)*exp(-2i*k*x);
A2_5_4_xx=diff( A2_5_4,x,2);
sum4_left=-simplify(4*k^2*A2_5_4+A2_5_4_xx);
param4_left=coeffs(expand(sum4_left/exp(-2i*k*x)),x);
[c4_1,c4_2,c4_3,c4_4,c4_5,c4_6,c4_7]=solve([param4_left(1)==param4_right(1),param4_left(2)==param4_right(2),param4_left(3)==param4_right(3),param4_left(4)==param4_right(4),param4_left(5)==param4_right(5),param4_left(6)==param4_right(6),param4_left(7)==param4_right(7)],[c4_1,c4_2,c4_3,c4_4,c4_5,c4_6,c4_7]);
A2_5_4=simplify(subs(A2_5_4));
%% 
% Eq5
sum5_right=(beta/2)*simplify(A1_4_c_x*A1_4_c_xx+A1con_4_c_x*A3_4_b_xx+A1con_4_c_xx*A3_4_b_x+A2con_4_b_x*A4_4_b_xx+A2con_4_b_xx*A4_4_b_x+A2con_4_c_x*A4_4_a_xx+A2con_4_c_xx*A4_4_a_x);
param5_right=coeffs(sum5_right/(exp(-2i*k*x)),x);
A2_5_5=(c5_1*x^9+c5_2*x^8+c5_3*x^7+c5_4*x^6+c5_5*x^5+c5_6*x^4+c5_7*x^3+c5_8*x^2+c5_9*x^1)*exp(-2i*k*x);
A2_5_5_xx=diff( A2_5_5,x,2);
sum5_left=-simplify(4*k^2*A2_5_5+A2_5_5_xx);
param5_left=coeffs(expand(sum5_left/exp(-2i*k*x)),x);
[c5_1,c5_2,c5_3,c5_4,c5_5,c5_6,c5_7,c5_8,c5_9]=solve([param5_left(1)==param5_right(1),param5_left(2)==param5_right(2),param5_left(3)==param5_right(3),param5_left(4)==param5_right(4),param5_left(5)==param5_right(5),param5_left(6)==param5_right(6),param5_left(7)==param5_right(7),param5_left(8)==param5_right(8),param5_left(9)==param5_right(9)],[c5_1,c5_2,c5_3,c5_4,c5_5,c5_6,c5_7,c5_8,c5_9]);
A2_5_5=simplify(subs(A2_5_5));
%% 
% Eq6
sum6_right=(beta/2)*simplify(A2con_4_c_x*A4_4_b_xx+A2con_4_c_xx*A4_4_b_x);
param6_right=coeffs(sum6_right/(exp(-2i*k*x)),x);
A2_5_6=(c6_1*x^11+c6_2*x^10+c6_3*x^9+c6_4*x^8+c6_5*x^7+c6_6*x^6+c6_7*x^5+c6_8*x^4+c6_9*x^3+c6_10*x^2+c6_11*x^1)*exp(-2i*k*x);
A2_5_6_xx=diff( A2_5_6,x,2);
sum6_left=-simplify(4*k^2*A2_5_6+A2_5_6_xx);
param6_left=coeffs(expand(sum6_left/exp(-2i*k*x)),x);
C6=solve([param6_left(1)==param6_right(1),param6_left(2)==param6_right(2),param6_left(3)==param6_right(3),param6_left(4)==param6_right(4),param6_left(5)==param6_right(5),param6_left(6)==param6_right(6),param6_left(7)==param6_right(7),param6_left(8)==param6_right(8),param6_left(9)==param6_right(9),param6_left(10)==param6_right(10),param6_left(11)==param6_right(11)],[c6_1,c6_2,c6_3,c6_4,c6_5,c6_6,c6_7,c6_8,c6_9,c6_10,c6_11]);
A2_5_6=simplify(subs(A2_5_6,C6));