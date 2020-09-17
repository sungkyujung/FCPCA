function cc = FCCCA(pc,lambda,k)
% FCCCA Functional Combined CCA. 
%
% cc = FCPCA(pc,lambda)
%      pc is the struct from FCPCA. This option skips the alignment
%      procedure by using the aligned result in pc. (For fast browing of
%      different options of c)
%      set lambda = 0.01 
%      set k = 3
% By Sungkyu Jung.
%
% See also FCPCA.m

if nargin < 1
    help FCCCA.m; return;
end

if nargin < 2
lambda = 0.01;
k = 3;
end


[d, n] = size(pc.data);

Sall = cov([pc.y ; pc.x]') + lambda * eye(2*d-1) ; 

S11 = Sall(1:d,1:d); 
S22 = Sall((d+1):end, (d+1):end);
S12 = Sall(1:d,(d+1):end); 

[V1,D1]=eig(S11); S11mh = ( V1 * diag(1./sqrt(diag(D1))) * V1' ) ; 
[V2,D2]=eig(S22); S22mh = ( V2 * diag(1./sqrt(diag(D2))) * V2' ) ; 
[u, dd, v]= svd(  S11mh  * S12 * S22mh, ...
     'econ');

cc.y = S11mh * u(:,1:k); 
cc.x = S22mh * v(:,1:k);
cc.canrho = diag(corr(pc.y' * cc.y, pc.x' * cc.x)); 

cc.scorey = pc.y' * cc.y ; 
cc.scorex = pc.x' * cc.x ;

end
 









