function  [P0,L0,Q0,X0,Z0,Y0,M0,Dc0,copt,m]=staticOutput_mincx(P0,L0,Q0,X0,Z0,Y0,M0,Dc0,A,Ad,B,C,Cd,dM,dm)
[numAi,numAj]=size(A);
[numAdi,numAdj]=size(Ad);
[numBi,numBj]=size(B);
[numCi,numCj]=size(C);
[numCdi,numCdj]=size(Cd);

setlmis([])
%定义矩阵变量
Dc=lmivar(2,[numBj,numCi]);
P=lmivar(1,[numAi,1]);
L=lmivar(1,[numAi,1]);
Q=lmivar(1,[numAi,1]);
X=lmivar(1,[numAi,1]);
Z=lmivar(1,[numAi,1]);
Y=lmivar(1,[numAi,1]);
M=lmivar(1,[numAi,1]);
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
m=decnbr(lmisys);
c=zeros(m,1);

    for j=1:m
        [Pj,Lj,Qj,Xj,Zj,Yj,Mj,Dcj]=defcx(lmisys,j,P,L,Q,X,Z,Y,M,Dc);
        c(j)=trace(Pj*L0+P0*Lj+Zj*M0+Z0*Mj);    %注意此时的参数P0,L0,Q0,X0,Z0,Y0,M0来自主程序，不是本函数生成，
    end
    
[copt,xopt]=mincx(lmisys,c);
P0=dec2mat(lmisys,xopt,P);
L0=dec2mat(lmisys,xopt,L);
Q0=dec2mat(lmisys,xopt,Q);
X0=dec2mat(lmisys,xopt,X);
Z0=dec2mat(lmisys,xopt,Z);
Y0=dec2mat(lmisys,xopt,Y);
M0=dec2mat(lmisys,xopt,M);
Dc0=dec2mat(lmisys,xopt,Dc);
%本函数的返回参数自动设置为第k+1步的参数
end