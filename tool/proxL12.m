function [x] = proxL12(y,alpha)
%x=proxL12(y,alpha)=argmin_x alpha(||x||_1-||x||_2)+1/2||x-y||_2^2
%   y>=0
%
dim=size(y);
if max(y)==0
    x=zeros(dim);
elseif  (0<max(y))&&(max(y)<=alpha)
      x=zeros(dim); 
     max_index = find(y == max(y));
     x(max_index(1))=1;
elseif max(y)>alpha
    z=max(y-alpha,0);
    x=(norm(z)+alpha)/(norm(z))*z;
end
end

