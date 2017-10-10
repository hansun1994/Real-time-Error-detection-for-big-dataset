clc
clear all
%load data1_y1_fault;
%load y1
%load y1_test_LOS
load ans

%% 选取对象
e=29;

for j=3:6
    y_tra=[];
    for u=1:j
        y_tra(u,:)=ans(1,e).CGMtrain(u:288-j+u,1);
    end

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
    
    t_main_stat(1,j)=min(t_stat);
end

%% 开始协整
    [C_stat,I_stat]=min(t_main_stat);
    y_tra=[];
    for u=1:I_stat
        y_tra(u,:)=ans(1,e).CGMtrain(u:288-I_stat+u,1);
    end

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

%% 选取胰岛素与饮食对象，并平滑
[y1_arx,data_injection_arx,data_CHO_arx]=ARX_pre(y_tra,ans(1,e).bolus(1:289-I_stat),ans(1,e).CHO(1:289-I_stat),1,1);
data_CHO_arx=data_CHO_arx';
data_injection_arx=data_injection_arx';

y_CHO_tra=data_CHO_arx;
y_injection_tra=data_injection_arx;

%% 检查平滑后的胰岛素和饮食数据的协整阶数
[CHO_Hvalue CHO_Delay CHO_Pvalue]=Judge_stationary(y_CHO_tra,Del);
[inject_Hvalue inject_Delay inject_Pvalue]=Judge_stationary(y_injection_tra,Del);

%% 将平滑后的胰岛素和饮食数据与t_main合并，构建新的变量
y_vec=[];
y_vec(1,:)=t_main;
y_vec(2,:)=y_CHO_tra;
y_vec(3,:)=y_injection_tra;

%% 验证是否平稳
[vec_Hvalue vec_Delay vec_Pvalue]=Judge_stationary(y_vec,Del);

%% 进行协整
p=3;              % The order of difference
deltay            = y_vec(:,p+1:end)-y_vec(:,p:end-1);      % delta y(k)
yk_1              = y_vec(:,p:end-1);                        % y(k-1)
Deltayk           = [];                                       % to store delta y(k-1) to delta y(k-p)
for i=1:p-1
    deltayk        = y_vec(:,p-i+1:end-i)-y_vec(:,p-i:end-1-i);
    Deltayk       = [Deltayk;deltayk];
end

[deltaHvalue deltaDelay deltaPvalue]=Judge_stationary(deltay,Del);

A                 = (deltay*Deltayk')*pinv(Deltayk*Deltayk'); 
R0                = deltay - A*Deltayk;

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
Eigenvevalue_vec      = [];
for s             = 1:size(SD,1)
    index_sum(s)  = sum(index(s:end));
    Eigenvevalue_vec  = [Eigenvevalue_vec V(:,order(s))];
end
t_Delay=[];
for s             = 1:size(SD,1)
    t_stat_vec(s)     = Judge_stationary1(Eigenvevalue_vec(:,s)'*y_vec,8);
end
[C_vec,I_vec]=min(t_stat_vec);
y_tra=Eigenvevalue_vec(:,s)'*y_vec;
[h, pvalue, stat, cValue, reg]= adftest(y_tra(1,:),'model','ARD','lags',5)

Eigenvevalue_vec(2,:)=-abs(Eigenvevalue_vec(2,:));

figure
t_main_vec=Eigenvevalue_vec(:,I_vec)'*y_vec;
%t_main_vec=Eigenvevalue_vec(:,1)'*y_vec;
plot(t_main_vec)
Del=8;
[tm_vec_Hvalue tm_vec_Delay tm_vec_Pvalue]=Judge_stationary(t_main_vec,Del);

%% 训练数据和测试数据的T2计算

y_vec_b=[];

[y1_arx,y_vec_b(3,:),y_vec_b(2,:)]=ARX_pre(y_tra,ans(1,e).bolus(289:577-I_stat),ans(1,e).CHO(289:577-I_stat),1,1);

%[p,cp,cpk]=capable(t_tra,[-3,3])
for i=1:11
    for u=1:I_stat
        y_test_error(u,:)=ans(1,e).CGMtest_bias(u:288-I_stat+u,i);
    end
    t_test_error_b=Eigenvevalue(:,I)'*y_test_error;
    y_vec_b(1,:)=t_test_error_b;
    t_vec_b=Eigenvevalue_vec(:,I_vec)'*y_vec_b;
    y_test=ans(1,e).CGMtest;
    [UUCL,max_rm_tra,rs_mean]=shewhart_test(t_main_vec',t_vec_b',y_test,ans(1,e).CGMtest_bias(:,i));
end
