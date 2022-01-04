% **************************
% 利用Hebb规则计算权值矩阵
% P=[p1,p2...,pn]，T=[t1,t2,..,tn]
% W=alpha*t1*p1'+alpha*t2*p2'+..+alpha*tn*pn'
% **************************
clc;clear;close all;
p1=[1;1;-1;1];p2=[1;-1;1;1];p3=[-1;-1;-1;1];
t1=p1;t2=p2;t3=p3;
P=[p1,p2,p3];                                                           %输入矩阵P[Q×N]
T=[t1,t2,t3];                                                             %期望输出矩阵T[S×N]
[S,~]=size(T);                                                       %期望矩阵的行数
[Q,~]=size(P);                                                      %输入矩阵的行数
Wold=zeros(S,Q);                                                %初始化权值矩阵
alpha=0.25;
for i=1:1:3
    a=purelin(Wold*P(:,i));
    Wnew=Wold+alpha*(T(:,i)-a)*P(:,i)';
    Wold=Wnew;
end
W=Wold;
hardlims(W*p1)
hardlims(W*p2)
hardlims(W*p3)