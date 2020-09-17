function pc = FCPCA(varargin)
% FCPCA Functional Combined PCA.
%
% Input :
% data    - d x n matrix of data, each column vector is
%                      a "d-dim digitized curve"
% grid    - a vector of length d, recording the grid points of the data
%           grid can be empty by setting [].
%
% pc = FCPCA(data) returns the combined principal components, where
%             data is the d x n matrix of data, each column vector is
%                      a "d-dim digitized curve"
%
% pc = FCPCA(data,grid)
%      grid is a vector of length d, recording the grid points of the data
%      grid can be empty by setting [].
%
% pc = FCPCA(data,grid,c) to supply additional parameters in estimation of
%       scaler c, where
%       c.nComp (default = 3) is the number of components `m' used in the
%       evaluation of the reconstruction errors, and
%       c.distanceoption  = 'L2';sum of L2-distance;
%                         = 'L1';sum of L1-distance;
%                         = 'FR'; sqrt of the squared sum of Fisher-Rao
%                         distance
%                         = 0;   default - sqrt of the squared sum of L2-distance
%
% pc = FCPCA(data,grid,c,k)
%      interpolate data on the domain of grid, but with k equi-spaced
%      points. ( k = 101 recommended for a balance of accuracy and speed)
%      set k = 0 to not interpolate.
%
% pc = FCPCA(pc,grid,c)
%      pc is the struct from FCPCA. This option skips the alignment
%      procedure by using the aligned result in pc. (For fast browing of
%      different options of c)
%
% By Sungwon Lee, Sungkyu Jung.
%
% See also FCPCAvis.m

if nargin < 1
    help FCPCA.m; return;
end
data = varargin{1};

% load SRVF FDA files (retrieved from http://ssamg.stat.fsu.edu/software
% 2015/08/08.
addpath(genpath('SRVF_FDA'));
 
if isstruct(data);
    pc = data;
    data = pc.data;
    [d_orig, n] = size(data);
    
    
    % set default parameters
    nComp = 3;
    dist = 0;
    
    if nargin == 2;
        gridpoints = varargin{2};
    if isempty(gridpoints) ; gridpoints = 1:d_orig ; end;
        creategrid = false; %k = 0;
    end
    
    if nargin >= 3;
        gridpoints = varargin{2};
    if isempty(gridpoints) ; gridpoints = 1:d_orig ; end;
        c = varargin{3};
        if isfield(c,'nComp'); nComp = c.nComp;end;
        if isfield(c,'distanceoption'); dist = c.distanceoption;end;
        creategrid = false; %k = 0;
     end
%     if nargin == 4;
%         grid = varargin{2};
%         c = varargin{3};
%         if isfield(c,'nComp'); nComp = c.nComp;end;
%         if isfield(c,'distanceoption'); dist = c.distanceoption;end;
%         k = varargin{4};
%         if k > 0 ; creategrid = true; gridsize = round(k); end;
%     end    
%    
%     if creategrid;
%         d = gridsize; % number of grid points
%         t = linspace(grid(1),grid(end),d);
%         f = zeros(d,n);
%         for i=1:n
%             f(:,i)=interp1(grid,data(:,i),t);
%         end
%     else
        t = 1:d_orig;
        f = data;
        d = d_orig;
%    end
    
    %% Transformation from gam to T_1 S.
    % step a: Compute SRVF of gam = phi
    h = 1/(d-1); % step size or the size of spacing
    phiS = zeros(d,n);
    for i = 1:n
        SRVF_func = sqrt(gradient(pc.gam(:,i),h));
        phiS(:,i) = SRVF_func / norm(SRVF_func);
    end
    % step b: Tangent approximation of SRVF on S
    muphi = geodmeanSm(phiS);
    rotmu = rotMat(muphi);
    x = LogNPd(rotmu*phiS);
    
    % centering y for easy computation.
    ycentered = bsxfun(@minus,pc.y,pc.mu_y);
    %% 2. Find empirical c, minimizing squared sum of L2 reconstruction errors.
    
    
    
    c = fminbnd(@(c) appxerror2(pc.data, ycentered, pc.mu_y, x, 1:nComp, rotmu, d , n, c, dist),0.001,500);
    
    
    %% 3. Return principal components (mean, score, and eigenstructure);
    
    g=[ycentered;c*x]; % g minus the mean
    
    %[u,score,latent]=pca(g');
    [u,s,v]=svd(g,'econ');
    %l = diag(s).^2 / (n-1);
    %score = s * v' / sqrt(n-1);
    
    % reconstruction of the mean:
    rec.phiS = (rotmu)' * ExpNPd(zeros(d-1,1)) ;
    rec.phiS2 = rec.phiS.*rec.phiS;
    rec.gam = [0 ;
        cumsum((  rec.phiS2(2:end-1,:) + rec.phiS2(1:end-2,:) ) /2 );
        1];
    cPCmean = interp1(linspace(0,1,d),pc.mu_y,rec.gam);
    
    %% 4. Return results
    %pc.data = data;
    %pc.grid = grid;
    %pc.mu_y = mu_y;
    %pc.mu_x = rotmu;
    pc.c =c;
    pc.latent = diag(s).^2 / (n-1);
    pc.score = s * v' ;
    pc.eigenframe = u;
    %pc.y = y;
    %pc.gam = gam;
    pc.mean = cPCmean;
    
else
    [d_orig,n]=size(data);
    
    % set default parameters
    nComp = 3;
    dist = 0;
    
    if nargin == 2;
        gridpoints = varargin{2};
        creategrid = false; %k = 0;
    end
    
    if nargin == 3;
        gridpoints = varargin{2};
        c = varargin{3};
        if isfield(c,'nComp'); nComp = c.nComp;end;
        if isfield(c,'distanceoption'); dist = c.distanceoption;end;
        creategrid = false; %k = 0;
    end
    if nargin == 4;
        gridpoints = varargin{2};
        c = varargin{3};
        if isfield(c,'nComp'); nComp = c.nComp;end;
        if isfield(c,'distanceoption'); dist = c.distanceoption;end;
        k = varargin{4};
        if k > 0 ; creategrid = true; gridsize = round(k); end;
    end
    
    if isempty(gridpoints) ; gridpoints = 1:d_orig ; end;
    
    if creategrid;
        d = gridsize; % number of grid points
        t = linspace(gridpoints(1),gridpoints(end),d);
        f = zeros(d,n);
        for i=1:n
            f(:,i)=interp1(gridpoints,data(:,i),t);
        end
    else
        t = gridpoints;
        f = data;
        d = d_orig;
    end
    
    %% 1. Decomposition into x and y variations
    %[fn,qn,q0,fmean,mqn,gam,psi,stats] = time_warping(f,t);
    % input: f and t
    % q0: the original SRVF function of the input f computed by
    % q0(:,i) = gradient(f(:,i), binsize)./sqrt(abs(gradient(f(:,i), binsize))+eps);
    % q0 is the data on which the Fisher-Rao alignment is performed.
    % After alignment, we have
    % fn: aligned functions ( = y )
    % qn: SRVF of fn
    % gam: warping functions (as in fn = f * gam )
    % psi: SRVF of gam (not needed)
    
    lambda = 0.1; % small lambda - more warping, larger lambda - less warping.
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 0;
    [y,~,~,~,~,gamINV,~,~] = time_warping(f,t,lambda,option);
    
    % retrieve warping functions as in f = y * gam, by inverting gamINV.
    gam = zeros(d,n);
    for i = 1:n
        gam(:,i) = invertGamma(gamINV(:,i));
    end
    
    %% Transformation from gam to T_1 S.
    % step a: Compute SRVF of gam = phi
    h = 1/(d-1); % step size or the size of spacing
    phiS = zeros(d,n);
    for i = 1:n
        SRVF_func = sqrt(gradient(gam(:,i),h));
        phiS(:,i) = SRVF_func / norm(SRVF_func);
    end
    % step b: Tangent approximation of SRVF on S
    muphi = geodmeanSm(phiS);
    rotmu = rotMat(muphi);
    x = LogNPd(rotmu*phiS);
    
    % centering y for easy computation.
    mu_y = mean(y,2);
    ycentered = bsxfun(@minus,y,mu_y);
    
    %% 2. Find empirical c, minimizing squared sum of L2 reconstruction errors.
    c = fminbnd(@(c) appxerror2(f, ycentered, mu_y, x, 1:nComp, rotmu, d , n, c, dist),0.001,500);
    
    
    %% 3. Return principal components (mean, score, and eigenstructure);
    
    g=[ycentered;c*x]; % g minus the mean
    
    %[u,score,latent]=pca(g');
    [u,s,v]=svd(g,'econ');
    %l = diag(s).^2 / (n-1);
    %score = s * v' / sqrt(n-1);
    
    % reconstruction of the mean:
    rec.phiS = (rotmu)' * ExpNPd(zeros(d-1,1)) ;
    rec.phiS2 = rec.phiS.*rec.phiS;
    rec.gam = [0 ;
        cumsum((  rec.phiS2(2:end-1,:) + rec.phiS2(1:end-2,:) ) /2 );
        1];
    cPCmean = interp1(linspace(0,1,d),mu_y,rec.gam);
    
    %% 4. Return results
    pc.data = data;       % input data
    pc.grid = gridpoints; % input grid
    pc.mu_y = mu_y;       % hat{mu}_y (mean of aligned functions)
    pc.mu_x = rotmu;      % mean of the time-warping functions (as a rotation matrix)
    pc.c =c;              % chosen value of the weight "c"
    pc.latent = diag(s).^2 / (n-1);
    pc.score = s * v' ;
    pc.eigenframe = u;
    pc.y = y;             % aligned functions: "y" (representing amplitudes)
    pc.gam = gam;         % time-warping functions: "\gamma"
    pc.x = x;             % phase functions: "x"
    pc.mean = cPCmean;    % FCPCA mean 
    
    
end

end




function score = appxerror2(f, ycentered, mu_y, x, PCid, rotmu, d , n, ccand, dist)

g=[ycentered;ccand*x]; % g minus the mean

%[u,score,latent]=pca(g');
[u,s,v]=svd(g,'econ');
%l = diag(s).^2 / (n-1);

% reconstruction
rec.g = u(:,PCid) * s(PCid,PCid) * v(:,PCid)';
rec.y = bsxfun(@plus,rec.g(1:d,:),mu_y);

rec.x = rec.g((d+1):end,:) / ccand;
rec.phiS = rotmu' * ExpNPd(rec.x) ;
rec.phiS2 = rec.phiS.*rec.phiS;
rec.gam = [zeros(1,n) ;
    cumsum((  rec.phiS2(2:end-1,:) + rec.phiS2(1:end-2,:) ) /2 );
    ones(1,n)];

rec.f = zeros(d,n);
t1 = linspace(0,1,d);
for j=1:n
    rec.f(:,j)=interp1(t1,rec.y(:,j),rec.gam(:,j));
end

switch dist
    case 'L2' % sum of L2-distance
        score =  sum(  sqrt( sum((f - rec.f).^2)) ) ;
        % sqrt(   sum(sum((f - rec.f).^2)) )
    case 'L1'  % sum of L1-distance
        score =   sum(sum(   abs(f - rec.f) ));
    case  'FR' % sqrt of the squared sum of Fisher-Rao distance
        
        h = 1/(d-1); % step size or the size of spacing
        q0 = zeros(d,n);
        qf = zeros(d,n);
        for j = 1:n;
           q0(:,j) =  sign(f(:,j)) .* sqrt(abs(gradient( f(:,j),h)));
           qf(:,j) =  sign(rec.f(:,j)) .* sqrt(abs(gradient( rec.f(:,j),h)));
        end 
        
        score =  norm( q0 - qf , 'fro'); 
    otherwise   % sqrt of the squared sum of L2-distance
        score = norm( f - rec.f , 'fro');
        
end


end
















