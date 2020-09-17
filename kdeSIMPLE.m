function [kdeh xgrid]=kdeSIMPLE(rowvec)
%rowvec = imat(1,:);
n = length(rowvec);
xgrid = linspace(min(rowvec),max(rowvec),201);
nbin = length(xgrid);
dsd = std(rowvec) ; diqr = iqr(rowvec) ;
a = min([dsd; (diqr / 1.34)]) ; %  Using Silverman's notation
h = .9 * a * n^(-1/5) ; % Silverman's eqn. 3.31 'rule of thumb'


kdeh = repmat((rowvec ./ h),nbin,1) - repmat((xgrid' ./ h),1,n) ;
          %  efficient way to divide all dif's by h
          %  variable name "kde" is used to avoid creating too many biggies
        kdeh = exp(-(kdeh .^2) / 2) ;
          %  exponential part of Gaussian density
        kdeh = sum(kdeh,2)' ;
          %  sum part of kde, and make result a column vector
        kdeh = kdeh / (n * h * sqrt(2 * pi)) ;
          %  normalize, and mult by Gaussain density constant

end