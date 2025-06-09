function [X] = prox_nulclear_L(Y,alpha)
%X=prox_matrix_nulclear_F(Y,alpha)=argmin_X alpha(||X||_*-||X||_F)+1/2||X-Y||_F^2
[U,S,V]=svd(Y,'econ');
s=diag(S);
z= proxL12(s,alpha);
X=U*diag(z)*V';
end

