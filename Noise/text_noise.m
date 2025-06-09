
%% 参数设置
length = 5000000;
beta = 2;

%% 生成噪声
 y = randlap(length,beta); % 生成拉普拉斯噪声
 v=var(y)-2*beta^2 ;%近似等于0
%% 概率密度计算

% 估计概率密度
[yy,x]=ksdensity(y);

% 计算概率密度理论值
xx = transpose(-5:1e-1:5);
miu = 0;
probablity1 = 1 / ( 2*beta ) * exp( -abs(xx-miu) / beta );

%% 作图
figure(1)
plot (x,yy,'bo'); % 做概率分布折线图
xlabel('noise');ylabel('frequency');
hold on;
%plot (xx,probablity1,'LineWidth',2);
plot (x,yy,'LineWidth',2);

%% 参数设置
length = 5000000;
G_std = 2;

%% 生成噪声
 y = G_std*randn(length,1); % 生成高斯噪声
%% 概率密度计算

% 估计概率密度
[yy,x]=ksdensity(y);

% 计算概率密度理论值
xx = transpose(-5:1e-1:5);
probablity2 = normpdf(xx,0,G_std);

%% 作图
figure(2)
xlabel('noise');ylabel('frequency');
plot (x,yy,'bo'); % 做概率分布折线图
hold on;
%plot (xx,probablity2,'LineWidth',2);
plot (x,yy,'LineWidth',2);


%% 参数设置
length = 5000000;
G_std = 2;beta = 2;

%% 生成噪声
 y = G_std*randn(length,1)+randlap(length,beta); % 生成混合噪声
%% 概率密度计算

% 估计概率密度
[yy,x]=ksdensity(y);

% 计算概率密度理论值
xx = transpose(-5:1e-1:5);
miu = 0;

%% 作图
figure(3)
xlabel('noise');ylabel('frequency');
plot (x,yy,'bo'); % 做概率分布折线图
hold on;
plot (x,yy,'LineWidth',2);
