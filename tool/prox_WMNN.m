function [X] = prox_WMNN(Y,W,rho)
    [U,S,V] = svd(Y,'econ');
     diagS = diag(S);
     rho1=rho*ones(size(diagS)).*W;    
     diagSS=max(diagS-rho1,0);
     X = U*diag(diagSS)*V';   
end
    
    
