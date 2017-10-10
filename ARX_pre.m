function [xgsmo,xismo,xmsmo]=ARX_pre(xg,xi,xm,transfoption,samplT)
%选择对数据是否进行平滑处理
%transfoption选择是否进行平滑处理
 %samplT,一般取5或1min
a=xg;
a1=xi;
a2=xm;
if transfoption==1
    h1=tf([1],[750 55 1]);
    h2=tf([1],[450 55 1]);
    Insulintr=c2d(h1,samplT,'tustin');
    Mealtr=c2d(h2,samplT,'tustin');
    [yout1,xout1]=lsim(Insulintr,a1);%计算Insulintr对a1的响应
    [yout2,xout2]=lsim(Mealtr,a2);
    a1=yout1;
    a2=yout2;
end
xgsmo=a;
xismo=a1;
xmsmo=a2;
