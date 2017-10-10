%休哈特控制图
function [UCL_rs,rm_tra,rm]=shewhart_test(X,Y,y_test,y_test_error)
[m,n]=size(Y);%求行数和列数
[m1,n1]=size(X);
data=0;
for i=1:m1
 data=data+X(i,1);
end
 data_mean=data/m1; 
 
 rs=0;
for j=1:m1-1
    rs=rs+abs(X(j+1,1)-X(j,1));
end
 rs_mean=rs/(m1-1);
 
for k=1:m-1
    rm(k,1)=abs(Y(k+1)-Y(k));
end

for l=1:m1-1
    rm_tra(l,1)=abs(X(l+1)-X(l));
end

 alpha=1.02*(max(rm_tra)/rs_mean);
 
 CL_x=data_mean;
 UCL_x=data_mean+0.36*rs_mean;
 LCL_x=data_mean-0.36*rs_mean;
 CL_rs=rs_mean;
 UCL_rs=alpha*rs_mean;
 LCL_rs=0.02*rs_mean;
 
 show=0.06*abs(y_test-y_test_error);

 figure
 plot(rm); hold on;
 %plot(rm,'LineWidth',6); hold on;
 plot(show,'g');hold on;
 for j=1:m
  plot(j,UCL_rs,'--r','LineWidth',6);hold on;
  plot(j,CL_rs,'--b');hold on;
  plot(j,LCL_rs,'--r');hold on;
 end
 