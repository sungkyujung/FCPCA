function  f = FCCCAscore2function(cc, pc , iCC, inputscore)


if nargin > 3
    
    nn = length(inputscore);
    d = length(pc.mu_y);
    rec.y = bsxfun(@plus,cc.y(:,iCC) * inputscore,pc.mu_y);
    rec.x = cc.x(:,iCC) * inputscore  / pc.c;
else
    [d,nn] = size(pc.data);
    rec.y = bsxfun(@plus,cc.y(:,iCC) * cc.scorey(:,iCC)',pc.mu_y);
    rec.x = cc.x(:,iCC) * cc.scorex(:,iCC)'  ; %/ pc.c;
    
end

rec.phiS = (pc.mu_x)' * ExpNPd(rec.x) ;
rec.phiS2 = rec.phiS.*rec.phiS;
rec.gam = [zeros(1,nn) ;
    cumsum((  rec.phiS2(2:end-1,:) + rec.phiS2(1:end-2,:) ) /2 );
    ones(1,nn)];
f = zeros(d,nn);
t1 = linspace(0,1,d);
for j=1:nn
    f(:,j)=interp1(t1,rec.y(:,j),rec.gam(:,j));
end
