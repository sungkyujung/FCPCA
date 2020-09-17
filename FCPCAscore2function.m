function  f = FCPCAscore2function(pc, inputscore, iPC)

    
    [nPC, nn]=size(inputscore);
    if nargin > 2 && nPC == 1; % Assume standardized scores 
            PCiscore = inputscore;  
    else 
         iPC = 1:nPC; 
         PCiscore = inputscore; 
    end
     
    d = length(pc.mu_y);
    rec.g = (pc.eigenframe(:,iPC) * PCiscore);
    rec.y = bsxfun(@plus,rec.g(1:d,:),pc.mu_y);    
    rec.x = rec.g((d+1):end,:) / pc.c;
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
    
    
    
    