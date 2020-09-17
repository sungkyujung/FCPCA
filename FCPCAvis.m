function  FCPCAvis(pc,visualizeType)
%
% FCPCAvis(pc) visually represent the combined variation of amplitude and
% pahse, contained in pc (output of FCPCA2.m).
%
% FCPCAvis(pc,visualizeType), choose visualizeType among 1, 2, ..., 4
%
%   1. cPC plot (curvdat style) - default.
%   2. cPC scores plot
%   3. Comparison to FPCA, Vertical PCA and Horizontal PCA (enabled with data input).
%   4. cPC plot (y,x, and combined)
%
% By Sungwon Lee, Sungkyu Jung.
%
% See also FCPCA.m


%% initialize
itype = 1;
if nargin>1
    itype = visualizeType;
end

[d, nn] = size(pc.y);

t1 = linspace(0,1,d);

f = pc.data;
t = pc.grid;
Titlestring1 = 'Raw data';




%%
clf;
[~,ind]=sort(pc.score(1,:));
cmat = rainbowcol(nn);
cmat(ind,:) = cmat;

switch itype
    case 1 % cPC plot (curvdat style) - default.
        
        subplot(5,3,1); hold on;
        for i = 1:nn
            plot(t,f(:,i),'Color',cmat(i,:));
        end
        xlim([t(1), t(end)]); title(Titlestring1);
        ylimits = get(gca,'ylim');
        
        subplot(5,3,2);
        plot(t,pc.mean,'Color','k','Linewidth', 2);
        xlim([t(1), t(end)]); ylim([ylimits(1) ylimits(2)]);
        title('FCPCA mean');
        
        subplot(5,3,3)
        per = pc.latent / sum(pc.latent);
        nPCshow = ceil(find(cumsum(per)> .95,1)/10)*10;
        bar(1:nPCshow,per(1:nPCshow),'LineStyle',':','FaceColor',[0 0.9 0.9],'EdgeColor',[0 0.7 0.7]); hold on;
        plot(1:nPCshow, cumsum([per(1:nPCshow)]),'-m'); hold off;
        xlim([.5 nPCshow+0.5]);
        title('Component variances');
        
        colors7 = colormap(winter(7));
        for iPC = 1:4;
            %  compute PC i projection / +- 3 PC and SCORES.
            
            subplot(5,3, 3*(iPC) + 2 )
            nm = 3;
            n = nm * 2 + 1;
            rec.PCiProjection= FCPCAscore2function(pc, sqrt(pc.latent(iPC))*(-nm:1:nm), iPC);
            hold on;
            for i = 1:n
                hp{i} = plot(t, rec.PCiProjection(:,i)  ,'Color',colors7(i,:));
            end
            xlim([t(1),t(end)]);
            ylim([ylimits(1) ylimits(2)]);
            set(hp{nm+1},'Color',[0,0,0],'LineStyle','-','LineWidth',2);
            title(['Center +- PC' num2str(iPC) ', ' num2str(per(iPC)*100,'%2.2f') '%']);
            
            subplot(5,3, 3*(iPC) + 1 )
            rec.PCiProjection= FCPCAscore2function(pc, pc.score(iPC,:), iPC);
            hold on;
            for i = 1:nn
                plot(t, rec.PCiProjection(:,i)  ,'Color', cmat(i,:));
            end
            xlim([t(1),t(end)]);
            ylim([ylimits(1) ylimits(2)]);
            title(['PC' num2str(iPC) ' Proj.' ]);
            
            subplot(5,3,  3*(iPC) + 3 )
            [kde, xgrid]=kdeSIMPLE(pc.score(iPC,:));
            heights = linspace(max(kde),2*max(kde),nn);
            plot(xgrid, kde); hold on;
            scatter(pc.score(iPC,:),heights,10,cmat);
            hold off; title(['PC' num2str(iPC) ' Score' ]);
            
        end
        
        
    case 2 %   2. cPC scores plot
        
        per = pc.latent / sum(pc.latent);
        paramstruct.npcadiradd = 0;
        paramstruct.irecenter = 0;
        paramstruct.icolor = cmat;
        paramstruct.labelcellstr = {['PC1, ' num2str(per(1)*100,'%2.2f') '%'],...
            ['PC2, ' num2str(per(2)*100,'%2.2f') '%'],...
            ['PC3, ' num2str(per(3)*100,'%2.2f') '%']};
        paramstruct.titlecellstr = {'','FCPCA Scores',''};
        
if exist('scatplotSM.m'); 
        scatplotSM(pc.score(1:3,:),eye(3),paramstruct)
else
      plotmatrix(pc.score(1:3,:)' );
      title([paramstruct.labelcellstr{1} ', ' paramstruct.labelcellstr{2} ', ' paramstruct.labelcellstr{3}])
end
        
        
        
    case 4 %% cPC plot;  +- 3 std.
        
        
        per = pc.latent / sum(pc.latent);
        nm = 3;
        n = nm * 2 + 1;
        for iPC = 1:3
            colors = rainbowcol(n);
            % reconstruction
            rec.g = kron(pc.eigenframe(:,iPC),  sqrt(pc.latent(iPC))*(-nm:1:nm));
            rec.y = bsxfun(@plus,rec.g(1:d,:),pc.mu_y);
            rec.x = rec.g((d+1):end,:) / pc.c;
            rec.phiS = (pc.mu_x)' * ExpNPd(rec.x) ;
            rec.phiS2 = rec.phiS.*rec.phiS;
            rec.gam = [zeros(1,n) ;
                cumsum((  rec.phiS2(2:end-1,:) + rec.phiS2(1:end-2,:) ) /2 );
                ones(1,n)];
            rec.f = zeros(d,n);
            t1 = linspace(0,1,d);
            for j=1:n
                rec.f(:,j)=interp1(t1,rec.y(:,j),rec.gam(:,j));
            end
            
            % column 1 : amplitudes
            subplot(3,3,1 + 3*(iPC-1));hold on;
            for i = 1:n
                hp{i} = plot(t, rec.y(:,i) ,'Color',colors(i,:));
            end
            xlim([t(1),t(end)]);
            set(hp{nm+1},'Color',[0,0,0],'LineStyle','-','LineWidth',2);
            title([num2str(iPC) 'th PC in Y']);
            
            % colum 2 : phases
            subplot(3,3,2 + 3*(iPC-1));hold on;
            for i = 1:n
                hp{i} = plot(t, t(1) + (t(end)-t(1))*rec.gam(:,i) ,'Color',colors(i,:));
            end
            xlim([t(1),t(end)]);ylim([t(1),t(end)]);
            set(hp{nm+1},'Color',[0,0,0],'LineStyle','-','LineWidth',2);
            title([num2str(iPC) 'th PC in X']);
            
            % colum 3 : combined.
            subplot(3,3,3 + 3*(iPC-1));hold on;
            for i = 1:n
                hp{i} = plot(t, rec.f(:,i)  ,'Color',colors(i,:));
            end
            xlim([t(1),t(end)]);
            set(hp{nm+1},'Color',[0,0,0],'LineStyle','-','LineWidth',2);
            title(['Reconsructed, ' num2str(per(iPC)*100,'%2.2f') '%']);
        end
        
        
    case 3 % comparison with FPCA
        
        
         subplot(3,3,1); hold on;
        for i = 1:nn
            plot(t,f(:,i),'Color',cmat(i,:));
        end
        xlim([t(1), t(end)]); title(Titlestring1);
        ylimits = get(gca,'ylim');
        
        subplot(3,3,2); hold on;
        for i = 1:nn
            plot(t,f(:,i),'Color',1-0.3*(1-cmat(i,:)));
        end
                h1 = plot(t,pc.mean,'Color','k','Linewidth', 1.5);
                h2 = plot(t,mean(f,2),':','Color',[0.2 0.2 0.5],'Linewidth', 1.5);
        xlim([t(1), t(end)]); ylim([ylimits(1) ylimits(2)]);
        legend([h1 h2],'FCPCA','Func. PCA','Location','BestOutside');
        title('Means');
        
        % 1. FPCA projections
          [u1, score1, latent1] = princomp(f');
          perFPCA = cumsum(latent1)/ sum(latent1); 
          for iPC = 1:2
          subplot(3,3,6+ iPC); 
          rec.FPC = bsxfun(@plus,u1(:,iPC) * (score1(:,iPC)'),  mean(f,2) ); 
          
            hold on;
            for i = 1:nn
                plot(t, rec.FPC(:,i)  ,'Color', cmat(i,:));
            end
            xlim([t(1),t(end)]);
            ylim([ylimits(1) ylimits(2)]);
            title(['Func. PC' num2str(iPC) ' Proj.' ]);
          end
          
          subplot(3,3,6+3);  
          scatter(score1(:,1),score1(:,2),10,cmat)
          xlabel(['PC1, ' num2str(perFPCA(1)*100,'%2.2f') '%'])
          ylabel(['PC2, ' num2str(perFPCA(2)*100,'%2.2f') '%'])
          title('Func. PCA Scores');
                  
        % 2. FCPCA projections
           for iPC = 1:2
          subplot(3,3,3+ iPC); 
            rec.PCiProjection= FCPCAscore2function(pc, pc.score(iPC,:), iPC);
            hold on;
            for i = 1:nn
                plot(t, rec.PCiProjection(:,i)  ,'Color', cmat(i,:));
            end
            xlim([t(1),t(end)]);
            ylim([ylimits(1) ylimits(2)]);
            title(['PC' num2str(iPC) ' Proj.' ]);
           end
           
        per = pc.latent / sum(pc.latent);
          subplot(3,3,3+3);  
          scatter(pc.score(1,:),pc.score(2,:),10,cmat)
          xlabel(['PC1, ' num2str(per(1)*100,'%2.2f') '%'])
          ylabel(['PC2, ' num2str(per(2)*100,'%2.2f') '%'])
          title('FCPCA Scores');
           
%            
%            % 3. Compute reconstruction errors; 
%            
%            rec.FPCAerror = zeros(10,1);
%            for ii = 1:10;
%                iPC = 1:ii;
%                rec.FPC =  bsxfun(@plus,u1(:,iPC) * (score1(:,iPC)'),  mean(f,2));
%                rec.FPCAerror(ii) = norm( f - rec.FPC  , 'fro');
%            end
%            
%            rec.FcPCAerror = zeros(10,1);
%            for ii = 1:10;
%                iPC = 1:ii;
%                rec.FcPC =  FCPCAscore2function(pc, pc.score(iPC,:),iPC);
%                rec.FcPCAerror(ii) = norm( f - rec.FcPC  , 'fro');
%            end
%             
%             subplot(3,3,3); hold on;
%             plot(1:10,rec.FcPCAerror,'Color','k','Linewidth', 1.5);
%             plot(1:10,rec.FPCAerror,':','Color',[0.2 0.2 0.5],'Linewidth', 1.5);
           
    otherwise; % do nothing
        
        
        
        
end










% 
% 
% 
% 
% %% cPC plot;  +- 3 std.
% %% PC scores plot;
% % color assignment by PC1
% figure;
% [~,ind]=sort(pc.score(1,:));
% nn = length(ind);
% cmat = rainbowcol(nn);
% cmat = cmat(ind,:);
% paramstruct.npcadiradd = 0;
% paramstruct.irecenter = 0;
% paramstruct.icolor = cmat;
% paramstruct.labelcellstr = {['PC1, ' num2str(per(1)*100,'%2.2f') '%'],...
%     ['PC2, ' num2str(per(2)*100,'%2.2f') '%'],...
%     ['PC3, ' num2str(per(3)*100,'%2.2f') '%']};
% paramstruct.titlecellstr = {'PC1','PC2','PC3'};
% scatplotSM(pc.score(1:3,:),eye(3),paramstruct)
% 
% 
% if nargin > 1;
%     h = figure;
%     
%     subplot(2,2,1); hold on;
%     for i = 1:nn
%         plot(t,f(:,i),'Color',cmat(i,:));
%     end
%     xlabel('t'); xlim([t(1), t(end)]); title('Raw data');
%     
%     subplot(2,2,2); hold on;
%     for i = 1:nn
%         plot(t,pc.y(:,i),'Color',cmat(i,:));
%     end
%     xlabel('t'); xlim([t(1), t(end)]);
%     title('Aligned data');
%     
%     subplot(2,2,3); hold on;
%     for i= 1:nn
%         plot(t,t(1) + (t(end)-t(1))*pc.gam(:,i),'Color',cmat(i,:));
%     end
%     xlabel('t');
%     xlim([t(1), t(end)]); ylim([t(1), t(end)]); title('warping functions');
% end
% hold off;
% 
% subplot(2,2,4); hold on;
% for i = 1:nn
%     plot(t,rec.f_orig(:,i),'Color',cmat(i,:));
% end
% xlabel('t'); xlim([t(1), t(end)]); title('Raw data');
% 
% 
% 
% 
% 
% 





