% **************************
% ����Hebb�������Ȩֵ����
% P=[p1,p2...,pn]��T=[t1,t2,..,tn]
% W=alpha*t1*p1'+alpha*t2*p2'+..+alpha*tn*pn'
% **************************
clc;clear;close all;
p1=[1;1;-1;1];p2=[1;-1;1;1];p3=[-1;-1;-1;1];
t1=p1;t2=p2;t3=p3;
P=[p1,p2,p3];                                                           %�������P[Q��N]
T=[t1,t2,t3];                                                             %�����������T[S��N]
[S,~]=size(T);                                                       %�������������
[Q,~]=size(P);                                                      %������������
Wold=zeros(S,Q);                                                %��ʼ��Ȩֵ����
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