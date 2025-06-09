
 addpath(genpath('RTC_GTNLN'));
 
 addpath(genpath('Noise'));
 addpath(genpath('data'));
 addpath(genpath('tool'));

 
 dataRoad = ['data/','guangzhou'];
 load(dataRoad);

 
disp('guangzhou_RTC_GTNLN');
traffictensor=guangzhou;
dim=size(traffictensor);

%The recovery results with various noisy patterns from noisy incomplete observations
missingway='R';
missingrate=0.5;
noisetensor=noise_tensor(dim,3,0);disp('50%LN-1+50%missing');
[~,~,~]=RTC_GTNLN(missingway,missingrate,noisetensor,traffictensor);

%The recovery results with various missing patterns from noisy incomplete observations
missingway='N';
missingrate=0.3;
noisetensor=noise_tensor(dim,3,0);disp('30%NM+70%noise');
[~,~,~]=RTC_GTNLN(missingway,missingrate,noisetensor,traffictensor);


 

