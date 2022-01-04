clear
clc
%矩阵参数
A=[0.9 0.5;0.8 0.1];Ad=[0.3 0;0.8 0.5];
B=[1;0.5];
C=[1 1;0 1];Cd=[1 0;1 1];
dM=11;
dm=3;
n=2;
N=100;%最大迭代次数
deta=1;
[numAi,numAj]=size(A);[numAdi,numAdj]=size(Ad);
[numBi,numBj]=size(B);
[numCi,numCj]=size(C);[numCdi,numCdj]=size(Cd);

setlmis([])
%定义矩阵变量
P=lmivar(1,[numAi,1]);
L=lmivar(1,[numAi,1]);
Q=lmivar(1,[numAi,1]);
X=lmivar(1,[numAi,1]);
Z=lmivar(1,[numAi,1]);
Y=lmivar(1,[numAi,1]);
M=lmivar(1,[numAi,1]);
Dc=lmivar(2,[numBj,numCi]);
%描述矩阵不等式
%(4)
lmiterm([-1 1 1 X],1,1);
lmiterm([-1 1 2 Y],1,1);
lmiterm([-1 2 2 Z],1,1);
%(20)
lmiterm([2 1 1 L],-1,1);
lmiterm([2 1 3 0],A);lmiterm([2 1 3 Dc],B,C);
lmiterm([2 1 4 0],Ad);lmiterm([2 1 4 Dc],B,Cd);
lmiterm([2 2 2 M],-1/(dM),1);
lmiterm([2 2 3 0],A);lmiterm([2 2 3 Dc],B,C);lmiterm([2 2 3 0],-1);
lmiterm([2 2 4 0],Ad);lmiterm([2 2 4 Dc],B,Cd);
lmiterm([2 3 3 P],-1,1);lmiterm([2 3 3 X],dM,1);lmiterm([2 3 3 Y],1,1,'s');lmiterm([2 3 3 Q],(dM-dm+1),1);
lmiterm([2 3 4 Y],-1,1);
lmiterm([2 4 4 Q],-1,1);
%(24)
lmiterm([-3 1 1 P],1,1);
lmiterm([-3 1 2 0],1);
lmiterm([-3 2 2 L],1,1);
lmiterm([-4 1 1 Z],1,1);
lmiterm([-4 1 2 0],1);
lmiterm([-4 2 2 M],1,1);

lmisys=getlmis;
matnbr(lmisys)   % 8个矩阵变量
lminbr(lmisys)   % 4个LMI矩阵不等式

%求上述不等式条件下的可行性解作为P0 L0 Q0 X0 Z0 Y0 M0
[tmin,xfeasp]=feasp(lmisys);
tmin
P0=dec2mat(lmisys,xfeasp,P);
L0=dec2mat(lmisys,xfeasp,L);
Q0=dec2mat(lmisys,xfeasp,Q);
X0=dec2mat(lmisys,xfeasp,X);
Z0=dec2mat(lmisys,xfeasp,Z);
Y0=dec2mat(lmisys,xfeasp,Y);
M0=dec2mat(lmisys,xfeasp,M);
Dc0=dec2mat(lmisys,xfeasp,Dc);
k=0;%迭代次数
%trace(P*Lk+Pk*L+Z*Mk+Zk*M)
%% 

while(true)
     [Pk,Lk,Qk,Xk,Zk,Yk,Mk,Dck,copt,m]=staticOutput_mincx(P0,L0,Q0,X0,Z0,Y0,M0,Dc0,A,Ad,B,C,Cd,dM,dm);
     ZZ=[-Pk^(-1) zeros(2,2) A+B*Dck*C Ad+B*Dck*Cd;
        zeros(2,2) (-1/dM)*Zk^(-1) A+B*Dck*C-eye(2) Ad+B*Dck*Cd;
        (A+B*Dck*C)' (A+B*Dck*C-eye(2))' -Pk+dM*Xk+Yk+Yk'+(dM-dm+1)*Qk -Yk;
        (Ad+B*Dck*Cd)' (Ad+B*Dck*Cd)' -Yk -Qk];
     t=abs((trace(Pk*Lk+Zk*Mk))-2*n)
%      t1=trace(Pk*Lk+Zk*Mk)
     [~,z]=eig(ZZ);
     z
     z_num=0;
    for j=1:1:length(z)
       if((z(j,j)>0)||(z(j,j)==0))
            break;
       end
      z_num=z_num+1; 
    end
     z_num
    if((z_num==length(z))&&t<deta)
        Dc0
        break;
     elseif(k<N)
        k=k+1
        P0=Pk;L0=Lk;Q0=Qk;X0=Xk;Z0=Zk;Y0=Yk;M0=Mk;Dc0=Dck;
    end
    if(k==N)
         Dc0
        break;
    end
end

