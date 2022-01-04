% 郑长波1，高杰. 应用线性矩阵不等式解决控制问题  沈阳师范大学学报(自然科学版) ，文献的matlab代码有很多错误
% 方法3：采用高会军文献中方法，主要区别在于更换判断是否有可行解的判断方法，就是判断整个Z矩阵的所有特征值是否都小于零
% 最后无解
clc;
clear all;

A=[1 0.25;0 1];
B=[0.0063;0.050];
C=[1 0];
D=0;

k=0


setlmis([])
P=lmivar(1,[2,1]);
Q=lmivar(1,[2,1]);
R=lmivar(1,[2,1]);
S=lmivar(1,[2,1]);
M=lmivar(1,[2,1]);
N=lmivar(1,[2,1]);
K=lmivar(2,[1,2]);
F=lmivar(2,[2,1]);

%(2)
lmi2=newlmi;
lmiterm([lmi2 1 1 S],1,1);
lmiterm([lmi2 1 1 P],-1,1);
lmiterm([lmi2 1 5 0],A);
lmiterm([lmi2 1 5 K],B,1);
lmiterm([lmi2 1 6 F],-1,C);
lmiterm([lmi2 2 2 S],-1,1);
lmiterm([lmi2 2 6 F],1,C);
lmiterm([lmi2 3 3 R],1,1);
lmiterm([lmi2 3 3 Q],-1,1);
lmiterm([lmi2 3 5 K],B,-1);
lmiterm([lmi2 3 6 F],1,C);
lmiterm([lmi2 3 6 0],-A);
lmiterm([lmi2 4 4 R],-1,1);
lmiterm([lmi2 5 5 M],-1,1);
lmiterm([lmi2 6 6 N],-1,1);

%(3)
lmi31=newlmi;
lmiterm([-lmi31 1 1 P],1,1);
lmiterm([-lmi31 1 2 0],1);
lmiterm([-lmi31 2 2 M],1,1);
lmi32=newlmi;
lmiterm([-lmi32 1 1 Q],1,1);
lmiterm([-lmi32 1 2 0],1);
lmiterm([-lmi32 2 2 N],1,1);

lmi_P=newlmi;
lmiterm([lmi_P 1 1 P],-1,1);
lmi_Q=newlmi;
lmiterm([lmi_Q 1 1 Q],-1,1);
lmi_R=newlmi;
lmiterm([lmi_R 1 1 R],-1,1);
lmi_S=newlmi;
lmiterm([lmi_S 1 1 S],-1,1);
lmi_M=newlmi;
lmiterm([lmi_M 1 1 M],-1,1);
lmi_N=newlmi;
lmiterm([lmi_N 1 1 N],-1,1);
lmisys=getlmis;

[tmin,xfeasp]=feasp(lmisys);
Pk=dec2mat(lmisys,xfeasp,P);
Qk=dec2mat(lmisys,xfeasp,Q);
Rk=dec2mat(lmisys,xfeasp,R);
Sk=dec2mat(lmisys,xfeasp,S);
Mk=dec2mat(lmisys,xfeasp,M);
Nk=dec2mat(lmisys,xfeasp,N);
Kk=dec2mat(lmisys,xfeasp,K);
Fk=dec2mat(lmisys,xfeasp,F);
n=decnbr(lmisys);
c=zeros(n,1);
 for j=1:n
        [Pj,Mj,Qj,Nj]=defcx(lmisys,j,P,M,Q,N);
        c(j)=trace(Pk*Mj+Pj*Mk+Qk*Nj+Qj*Nk);   %注意此函数的参数Pk,Qk,Mk,Nk来自上述feasp可行解的参数
 end