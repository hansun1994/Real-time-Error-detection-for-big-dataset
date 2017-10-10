clc
clear all
%load data1_y1_fault;
load y1
load y1_test_LOS
load all_data

%% 选取对象
e=1;

%y5=CGMfail(e,1).a;
y2=y1(e,:);
%y3=y2+sqrt(2)*randn(1,length(y2));
%y3=randn(1,length(y2));
y=y2(1:5:end);
y2_test=y1_test_LOS(e,:);
y3_test=y2_test(1:5:end);

n=5;
y_vector=[];
for i=1:n
    y_vector(i,:)=y(i:end-n+i);
    y_vector_error(i,:)=y3_test(i:end-n+i);
end

N_test=288;
y_tra=y_vector(:,1:N_test-n+1);
%y_tra=CGMfail(e,1).a';

%% 判断每个变量是否是平稳变量
Del=8;
[Hvalue Delay Pvalue]=Judge_stationary(y_tra,Del);

%% To calculate the eigenvectors and eigenvalues ()
p=3;              % The order of difference
deltay            = y_tra(:,p+1:end)-y_tra(:,p:end-1);      % delta y(k)
yk_1              = y_tra(:,p:end-1);                        % y(k-1)
Deltayk           = [];                                       % to store delta y(k-1) to delta y(k-p)
for i=1:p-1
    deltayk        = y_tra(:,p-i+1:end-i)-y_tra(:,p-i:end-1-i);
    Deltayk       = [Deltayk;deltayk];
end

[deltaHvalue deltaDelay deltaPvalue]=Judge_stationary(deltay,Del);

%% Calculate the residual sequence R0
A                 = (deltay*Deltayk')*pinv(Deltayk*Deltayk'); 
R0                = deltay - A*Deltayk;
%% Calculate the residual sequence R1
B                 = (yk_1*Deltayk')*pinv(Deltayk*Deltayk');
R1                = (yk_1-B*Deltayk);
S00               = R0*R0';
S01               = R0*R1';
S10               = R1*R0';
S11               = R1*R1';

[V,D]             = eig(S10*inv(S00)*S01,S11);
[SD,order]        = sort(diag(D),'descend');
for s=1:size(SD,1)
    index(s)      = -size(R0,2)*log(1-SD(s));
end
Eigenvevalue      = [];
for s             = 1:size(SD,1)
    index_sum(s)  = sum(index(s:end));
    Eigenvevalue  = [Eigenvevalue V(:,order(s))];
end
t_Delay=[];
for s             = 1:size(SD,1)
    t_stat(s)     = Judge_stationary1(Eigenvevalue(:,s)'*y_tra,8);
end
[C,I]=min(t_stat);
figure
t_main=Eigenvevalue(:,I)'*y_tra;
plot(t_main)
Del=8;
[tm_Hvalue tm_Delay tm_Pvalue]=Judge_stationary(t_main,Del);

%% 训练数据和测试数据的T2计算
y_tra=y_vector(:,1:N_test-n+1);
y_test=y_vector(:,N_test-n:end);
y_test_error=y_vector_error(:,N_test-n:end);

figure
subplot(311)
t_tra=Eigenvevalue(:,I)'*y_tra;
plot(t_tra)
subplot(312)
t_test=Eigenvevalue(:,I)'*y_test;
plot(t_test)
subplot(313)
t_test_error=Eigenvevalue(:,I)'*y_test_error;
plot(t_test_error)

[p,cp,cpk]=capable(t_tra,[-3,3])
[UUCL,max_rm_tra,rs_mean]=shewhart_test(t_tra',t_test_error');
