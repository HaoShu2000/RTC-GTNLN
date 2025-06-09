
function [X,MAE,RMSE] = RTC_GTNLN(missingway,missingrate,noisetensor,traffictensor)
%% default paremeters setting
X0=traffictensor;
Y=traffictensor+noisetensor;
dim = size(X0);
switch missingway
      case 'R'%%各个维度随机丢失 Random_Missing
            Pomega=round(rand(dim(1),dim(2),dim(3)) + 0.5 - missingrate);
      case 'N'   %%在某些时间点所有传感器全黑（Non_Random_Missing:tuba1)
            A=round(rand(dim(3), dim(2)) + 0.5 - missingrate);
            B= kron(A,ones(dim(1),1));
            Pomega=reshape(B,[dim(1),dim(2),dim(3)]);
end
Pomegac=1-Pomega;
tol        = 1e-3; 
max_iter   = 200;
max_mu     = 1e10;
detail     = 1;
rho        = 1.1;
mu         = 1e-6;
alpha=[1/3,1/3,1/3];
lambda =1/sqrt(max(dim(1),dim(2))*dim(3));
%% variables initialization
X            = Pomega.*Y;

D    =circvet2mat(num2vetT(1,dim(2))); 

G    = nmodeproduct(ones(dim),D,2); 

E    = zeros(dim);

K            = zeros(dim);

M      = zeros(dim); 

N      = zeros(dim);

Z{1}      = zeros(dim(1),dim(2)*dim(3)); 
Z{2}      = zeros(dim(2),dim(1)*dim(3)); 
Z{3}      = zeros(dim(3),dim(1)*dim(2)); 

Q{1}      = zeros(dim(1),dim(2)*dim(3)); 
Q{2}      = zeros(dim(2),dim(1)*dim(3)); 
Q{3}      = zeros(dim(3),dim(1)*dim(2)); 

[V,S]=eig(D'*D);

T=nmodeproduct(ones(dim),S,2);

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Ek = K;
    %% Update X
 
 
    MAE=(1/(prod(dim)))*sum(abs(X0-X),'all');
    RMSE=1/sqrt(prod(dim))*sqrt(sum((X0-X).^2,'all'));
    X= Pomega.*Y-K-E+N/mu+nmodeproduct(G-M/mu,D',2);
    
    X = nmodeproduct(X,V',2); 

    X=(X./(1+T));
  
    X = nmodeproduct(X,V,2); 


    %% Updata G
    H=nmodeproduct(X,D,2)+M/mu;
    for i=1:3
    H=H+Fold(Z{i}+Q{i}/mu,dim,i);
    end
    G=H/4;
 


    %% Updata Z_i
    
 for i=1:3
 
        [Z{i}] = prox_nulclear_L(Unfold(G,dim,i)-Q{i}/mu,alpha(i)/mu); 

 end
    %% Updata E
    E= prox_l1( Pomega.*Y-K-X+N/mu,lambda/mu);

      
  

 
    %% Update K 
    K          = Pomega.*Y-X-E+N/mu;
    K          = Pomegac.*K;
    
    %% Stop criterion
    dY   =  Pomega.*Y-X-E-K;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-K(:)));
    chg  = max([chgX chgE max(abs(dY(:)))]);
    if chg < tol
         break;
    end 
  
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 30) == 0
            err = norm(dY(:),'fro');
            disp(['iter= ' num2str(iter) ', mu=' num2str(mu) ...
                   ', chg=' num2str(chg) ...
                     ', err=' num2str(err) ',  MAE=' num2str(MAE) ...
                     ',  RMSE=' num2str( RMSE)]); 
        end
    end   
    
    %% Update mulipliers: Lambda, Gamma, and mu
  
    M= M+mu*(nmodeproduct(X,D,2)-G);
  
    N = N+mu*dY;
    
    for i=1:3
    Q{i} =Q{i} +mu*(Z{i}-Unfold(G,dim,i));
    end
        
    mu = min(rho*mu,max_mu);
end
end