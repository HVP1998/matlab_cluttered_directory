%***********************
% class1��p1=[1;1];p2=[1;2];
% class2��p3=[2;-1];p4=[2;0];
% class3��p5=[-1;2];p6=[-2;1];
% class4��p7=[-1;-1];p8=[-2;-2];
% ����class1��class2������t=0��class3��class4��������t=1
% Ȼ��class1��class3������t=0��class2��class4��������t=1
%***********************
t1=0;t2=1;                                                           %����t
W1old=[1,1];b1old=1;                                          %���ó�ʼW1��b1
W2old=[1,1];b2old=1;                                          %���ó�ʼW2��b2
p1=[1;1];p2=[1;2];
p3=[2;-1];p4=[2;0];
p5=[-1;2];p6=[-2;1];
p7=[-1;-1];p8=[-2;-2];
P={p1,p2,p3,p4,p5,p6,p7,p8};                               %����������������Ԫ��������
e=ones(1,8);
while(1)
    for i=1:1:8
        a=hardlim(W1old*cell2mat(P(i))+b1old);      %���ù�ʽ��a
        if(i<=4)
            e(i)=t1-a;
        else
            e(i)=t2-a;
        end
        W1new=W1old+e(i)*(cell2mat(P(i)))';           %�����µ�W����  
        W1old=W1new;                                          %����W����
        b1new=b1old+e(i);                                      %�����µ�b
        b1old=b1new;                                             %����b
    end
    if(e==zeros(1,8))                                              %��e�������������Ƴ�ѭ��
        e
        break;
    end
    e
end
while(1)
    for i=1:1:8
        a=hardlim(W2old*cell2mat(P(i))+b2old);
        if((i==1)||(i==2)||(i==5)||(i==6))
            e(i)=t1-a;
        else
            e(i)=t2-a;
        end
        W2new=W2old+e(i)*(cell2mat(P(i)))';
        W2old=W2new;
        b2new=b2old+e(i);
        b2old=b2new;
    end
    if(e==zeros(1,8))
        e
        break;
    end
    e
end
