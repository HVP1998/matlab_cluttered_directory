% **************************
% ���������½��������Ȩֵ����
% ���ڴ˴������ĸ�����������û�г˻������Կ��Խ��������ַֿ�����Сֵ
% **************************
clc;clear;close all;
x_a=[0,0,0,0];                                                       %��ʼ��
x_b=[2,3,-4,2];                                                      %��ֹ��
x_c=ones(1,4);                                                     %�ƽ�������   Ԥ�ȷ���ռ�
x_d=ones(1,4);                                                     %�ƽ�������   Ԥ�ȷ���ռ�
result=ones(1,4);                                                  %���ս������   Ԥ�ȷ���ռ�
e=1e-6;                                                                %��������
for i=1:1:4                                                            %�Ĵ�ѭ���ֱ��x1,x2,x3,x4��Сֵ
    if(x_a(i)<x_b(i))                                                  %�ж�a��b�Ǹ���Ϊ��Сֵ�ĸ���Ϊ��Сֵ
            min=x_a(i);
            max=x_b(i);
        else
            min=x_b(i);
            max=x_a(i);
    end
    while(1)
        x_c(i)=min+0.382*(max-min);                       %����ƽ�ָ����
        x_d(i)=min+0.618*(max-min);                       %����ƽ�ָ����
        switch(i)                                                       %i=nʱ����xn��һ��ĺ���ֵ
            case 1
                f_x_c=3*x_c(i)^2;
                f_x_d=3*x_d(i)^2;
            case 2
                f_x_c=(x_c(i)-1)^2;
                f_x_d=(x_d(i)-1)^2;
            case 3
                f_x_c=2*(x_c(i)+2)^2;
                f_x_d=2*(x_d(i)+2)^2;
            case 4
                f_x_c=x_c(i)^2;
                f_x_d=x_d(i)^2;
        end
%          f_x_c
%          f_x_d
        if( f_x_c>=f_x_d)                                          %���f_x_c>=f_x_d����ȥmin-c�Σ����f_x_c<f_x_d����ȥd-max��
            min=x_c(i);                                              %��ƽ���������Ϊȡֵ��Χ����
        else
            max=x_d(i);                                             %��ƽ���������Ϊȡֵ��Χ����
        end
        if(abs(min-max)<e)                                      %������㾫�ȵ������˳�
            break;
        end
    end
    x_result(i)=min;
end
 x_result
 f_result=3*x_result(i)^2+(x_result(i)-1)^2+2*(x_result(i)+2)^2+x_result(i)^2+10