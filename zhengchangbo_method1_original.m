% ֣����1���߽�. Ӧ�����Ծ��󲻵�ʽ�����������  ����ʦ����ѧѧ��(��Ȼ��ѧ��) �����׵�matlab�����кܶ����
% ����1��ѧ�����˻��������ڴ˻����ϸ�����Щ��ĵط�������trace���ʽ��Z��ת�����⣬�õķ�������ƪ���׵ģ�����޽⡣
clc;
clear all;

A=[1 0.25;0 1];
B=[0.0063;0.050];
C=[1 0];
D=0;

setlmis([])
P=lmivar(1,[2,1]);
Q=lmivar(1,[2,1]);
R=lmivar(1,[2,1]);
S=lmivar(1,[2,1]);
M=lmivar(1,[2,1]);
N=lmivar(1,[2,1]);
K=lmivar(2,[1,2]);
F=lmivar(2,[2,1]);

lmiterm([1 1 1 S],1,1);
lmiterm([1 1 1 P],-1,1);
lmiterm([1 1 5 0],A);
lmiterm([1 1 5 K],B,1);
lmiterm([1 1 6 F],-1,C);
lmiterm([1 2 2 S],-1,1);
lmiterm([1 2 6 F],1,C);
lmiterm([1 3 3 R],1,1);
lmiterm([1 3 3 Q],-1,1);
lmiterm([1 3 5 K],B,-1);
lmiterm([1 3 6 F],1,C);
lmiterm([1 3 6 0],-A);
lmiterm([1 4 4 R],-1,1);
lmiterm([1 5 5 M],-1,1);
lmiterm([1 6 6 N],-1,1);

lmiterm([-2 1 1 P],1,1);
lmiterm([-2 1 2 0],1);
lmiterm([-2 2 2 M],1,1);

lmiterm([-3 1 1 Q],1,1);
lmiterm([-3 1 2 0],1);
lmiterm([-3 2 2 N],1,1);

lmiterm([-4 1 1 P],1,1);
lmiterm([-5 1 1 Q],1,1);
lmiterm([-6 1 1 R],1,1);
lmiterm([-7 1 1 S],1,1);
lmiterm([-8 1 1 M],1,1);
lmiterm([-9 1 1 N],1,1);
lmisys=getlmis;

[tmin,xfeasp]=feasp(lmisys);
PP=dec2mat(lmisys,xfeasp,P);
QQ=dec2mat(lmisys,xfeasp,Q);
RR=dec2mat(lmisys,xfeasp,R);
SS=dec2mat(lmisys,xfeasp,S);
MM=dec2mat(lmisys,xfeasp,M);
NN=dec2mat(lmisys,xfeasp,N);
KK=dec2mat(lmisys,xfeasp,K);
FF=dec2mat(lmisys,xfeasp,F);
%% 
for i=1:100                                               %���õ���������Ϊ100
    n=decnbr(lmisys);
    c=zeros(n,1);
    
    for j=1:n
        [Pj,Mj,Qj,Nj]=defcx(lmisys,j,P,M,Q,N);
        c(j)=trace(PP*Mj+Pj*MM+QQ*Nj+Qj*NN); %ע��˺����Ĳ���PP,QQ,MM,NN��������feasp���н�Ĳ���
    end
    [copt,xopt]=mincx(lmisys,c);
    PPP=dec2mat(lmisys,xopt,P)
    QQQ=dec2mat(lmisys,xopt,Q)
    RRR=dec2mat(lmisys,xopt,R)
    SSS=dec2mat(lmisys,xopt,S)
    MMM=dec2mat(lmisys,xopt,M)
    NNN=dec2mat(lmisys,xopt,N)
    KKK=dec2mat(lmisys,xopt,K)
    FFF=dec2mat(lmisys,xopt,F)
    Z=[SSS-PPP,zeros(2),zeros(2),zeros(2),(A+B*KKK),-(FFF*C);
       zeros(2),-SSS,zeros(2),zeros(2),zeros(2),(FFF*C);
       zeros(2),zeros(2),RRR-QQQ,zeros(2),-(B*KKK),FFF*C-A;
       zeros(2),zeros(2),zeros(2),-RRR,zeros(2),zeros(2);
       (A+B*KKK)',zeros(2),(-B*KKK)',zeros(2),-inv(PPP),zeros(2);
       (-FFF*C)',(FFF*C)',(FFF*C-A)',zeros(2),zeros(2),-inv(QQQ)];
   Y=eig(Z)   %��Ϊ����Z�����ά��Ϊ12��12��Z����������ֵ��Y������һ����12����YΪ12��1������
   length(Y);  %length(Y)=12

   i2=0;
   for i1=1:length(Y)   % i1��ȡֵ��Χ��1��12
       if(Y(i1,1)<0)    %Y(i1,1)��Y�ĵ�i1������ֵ,�жϵ�i1������ֵ�Ƿ�С��0
           i2=i2+1;
       end
   end
   if(i2==length(Y))
       break;
   end
   i
   PP=PPP
   QQ=QQQ;
   RR=RRR;
   SS=SSS;
   MM=MMM;
   NN=NNN;
   KK=KKK;
   FF=FFF
   copt
end
if(i==100)
    disp('There is no result');
end

%{
��100�ε����ж�Z��12������ֵ
Y =
  -12.6568
   -7.9144
   -4.3370
   -2.4359
   -0.2801
   -0.1367
   -0.0386
   -0.0092
   -0.0086
   -0.0081
   -0.0011
    0.0092
����ֵ�����1��������������Z�����У��Ѿ��ﵽ����������������ѭ�����޽⡣
ע��˷�����Method3��100����ʱ�õ���12������ֵ��ͬ��˼��Ϊ�β�ͬ������ϰ5�лش��������
%}