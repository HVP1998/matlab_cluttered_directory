% **************************
% 利用最速下降法求计算权值矩阵
% 由于此处函数的各个变量都是没有乘积的所以可以将各个部分分开算最小值
% **************************
clc;clear;close all;
x_a=[0,0,0,0];                                                       %初始点
x_b=[2,3,-4,2];                                                      %终止点
x_c=ones(1,4);                                                     %黄金分配点左   预先分配空间
x_d=ones(1,4);                                                     %黄金分配点右   预先分配空间
result=ones(1,4);                                                  %最终结果向量   预先分配空间
e=1e-6;                                                                %收敛精度
for i=1:1:4                                                            %四次循环分别对x1,x2,x3,x4求极小值
    if(x_a(i)<x_b(i))                                                  %判断a，b那个作为最小值哪个作为最小值
            min=x_a(i);
            max=x_b(i);
        else
            min=x_b(i);
            max=x_a(i);
    end
    while(1)
        x_c(i)=min+0.382*(max-min);                       %计算黄金分割点左
        x_d(i)=min+0.618*(max-min);                       %计算黄金分割点右
        switch(i)                                                       %i=n时计算xn那一项的函数值
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
        if( f_x_c>=f_x_d)                                          %如果f_x_c>=f_x_d则舍去min-c段，如果f_x_c<f_x_d则舍去d-max段
            min=x_c(i);                                              %令黄金分配点左作为取值范围下限
        else
            max=x_d(i);                                             %令黄金分配点右作为取值范围上限
        end
        if(abs(min-max)<e)                                      %如果计算精度到达则退出
            break;
        end
    end
    x_result(i)=min;
end
 x_result
 f_result=3*x_result(i)^2+(x_result(i)-1)^2+2*(x_result(i)+2)^2+x_result(i)^2+10