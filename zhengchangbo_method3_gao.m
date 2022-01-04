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
%注意k=0,所以此时Pk即为P0，也就是第一个可行解,大家可以看一下，此时是可以找到一组可行解P0,Q0,R0,S0,M0,N0,K0,F0
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 算法
n=decnbr(lmisys);
c=zeros(n,1);

iterate_no=101; %设置迭代最大次数为100
ttt=1;
while ttt>0  %ttt先设置为1，后面取为特征值最大值作为判断条件
    k=k+1
    if k>=iterate_no  
        break
    else
    for j=1:n
        [Pj,Mj,Qj,Nj]=defcx(lmisys,j,P,M,Q,N);
        c(j)=trace(Pk*Mj+Pj*Mk+Qk*Nj+Qj*Nk);   %注意此函数的参数Pk,Qk,Mk,Nk来自上述feasp可行解的参数
    end
[copt,xopt]=mincx(lmisys,c)
Popt=dec2mat(lmisys,xopt,P)
Qopt=dec2mat(lmisys,xopt,Q)
Ropt=dec2mat(lmisys,xopt,R)
Sopt=dec2mat(lmisys,xopt,S)
Mopt=dec2mat(lmisys,xopt,M)
Nopt=dec2mat(lmisys,xopt,N)
Kopt=dec2mat(lmisys,xopt,K)
Fopt=dec2mat(lmisys,xopt,F)
    
Pk=Popt;
Qk=Qopt;
Rk=Ropt;
Sk=Sopt;
Mk=Mopt;
Nk=Nopt;
Kk=Kopt;
Fk=Fopt;

Z=[Sk-Pk,zeros(2),zeros(2),zeros(2),(A+B*Kk),-(Fk*C);
   zeros(2),-Sk,zeros(2),zeros(2),zeros(2),(Fk*C);
   zeros(2),zeros(2),Rk-Qk,zeros(2),-(B*Kk),Fk*C-A;
   zeros(2),zeros(2),zeros(2),-Rk,zeros(2),zeros(2);
   (A+B*Kk)',zeros(2),(-B*Kk)',zeros(2),-inv(Pk),zeros(2);
   (-Fk*C)',(Fk*C)',(Fk*C-A)',zeros(2),zeros(2),-inv(Qk)];
Y=eig(Z)           %12个特征值
ttt=max(eig(Z))  %特征值中的最大值
end
end

k
PP=Pk
QQ=Qk;
RR=Rk;
SS=Sk;
MM=Mk;
NN=Nk;
KK=Kk;
FF=Fk
copt

if(k==iterate_no)
    disp('There is no result');
end

%{
大家看一下，
第1次迭代判断Z的12个特征值时
Y =
   -7.8194
   -6.5071
   -2.1854
   -1.7920
   -0.1675
   -0.0446
   -0.0243
   -0.0102
   -0.0045
   -0.0043
    0.0458
    0.1185
ttt =
    0.1185
特征值的最后两个是正数，所以Z不可行，同时最大特征值ttt大于0，继续循环，使得都要小于0

第100次迭代判断Z的12个特征值时
Y =
 -106.3697
  -79.8424
  -13.2215
   -2.4815
   -1.7870
   -0.7255
   -0.3389
   -0.1735
   -0.1095
   -0.0143
   -0.0007
    0.0026
ttt =
    0.0026
特征值的最后1个是正数，所以Z不可行，同时最大特征值ttt仍旧大于0，已经达到最大迭代次数，不再循环，无解。
注意此方法和Method1第100迭代时得到的12个特征值不同，思考为何不同？在练习5中回答这个问题

此外，可以验证当最大迭代次数设置为1000时，最后一个特征值为2.5701e-004，无限接近0，但还是大于0,仍然无解。
%}