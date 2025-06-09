

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
probablity = 1 / ( 2*beta ) * exp( -abs(xx-miu) / beta );

%% 作图
figure;xlabel('x');ylabel('PDF');
plot (x,yy,'bo'); % 做概率分布折线图
hold on;
plot (xx,probablity,'LineWidth',2);
legend('PDF');

