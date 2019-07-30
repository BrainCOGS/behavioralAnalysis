function stats = clusterByOptoEffect(lsrPerf,cfg,plotFlag)

% stats = clusterByOptoEffect(lsrPerf,cfg,plotFlag)
% clusters locations based on laser effect
% lsrPerf is output structure of analyzeConcatLsrLog()
% cfg is optional analysis parameter structure
% plotFlag is boolean to plot data

if nargin < 2
  %% analysis config
  cfg.clustWhat        = {'percCorrect','bias_abs','speed','excessTravel','unusualEvents'};
  cfg.clustMaxPC       = 3; % number of PCs to use for clustering
  cfg.clustMaxK        = 5; % maximal number of clusters to test
  cfg.distMeasure      = 'eucl'; % 'eucl' or 'correlation'
  cfg.fixedK           = false; % fixed number of clusters? (= maxK)
  %%
end
if nargin < 3; plotFlag = false; end
stats.cfg = cfg;

%% collect desired effect sizes for clustering and calculate correlations
grid = lsrPerf{1}.info.grid;
ML   = flipud(cellfun(@(x)(x(2,1)),grid)'); % flip for convenience to have visual areas first
AP   = flipud(cellfun(@(x)(x(2,2)),grid)');

effectSize = zeros(numel(ML),numel(cfg.clustWhat));
for iVec = 1:numel(cfg.clustWhat)
  pvals      = retrievePvals(lsrPerf,cfg.clustWhat{iVec},'boot');
  effectSize(:,iVec) = flipud(pvals.effectSize(1:2:end,:));
end

cc = corr(effectSize');%effectSize';%

%% do PCA and clustering
[pcload,~,ev]    = pca(cc);
stats.fracVarPCs = sum(ev(1:cfg.clustMaxPC))/sum(ev);
[stats.clustID,stats.distMat,stats.lnk,stats.CHindex] ...
                 = hierClustering(pcload(:,1:cfg.clustMaxPC),cfg.distMeasure,cfg.clustMaxK,cfg.fixedK);
stats.nK         = numel(unique(stats.clustID));
clustID          = stats.clustID;
stats.dataMat    = effectSize;
stats.corrMat    = cc;
stats.ML         = ML;
stats.AP         = AP;
stats.fracTopPCs = sum(ev(1:cfg.clustMaxPC)./sum(ev))*100;

%% do stats to compare significance of effect differences between clusters
% one-way ANOVA separately for each weight group
stats.labels     = cellfun(@removeUnderscores,cfg.clustWhat,'uniformoutput',false);

for iPred = 1:numel(stats.labels)
  datavec = [];
  groupid = [];
  for iK = 1:stats.nK
    datavec = [datavec; stats.dataMat(stats.clustID == iK,iPred)];
    groupid = [groupid; ones(sum(stats.clustID == iK),1).*iK];
  end
  [stats.ANOVA(iPred).pval,stats.ANOVA(iPred).table,stats.ANOVA(iPred).anovaStats] ...
                          = anovan(datavec,groupid,'display','off');
  stats.ANOVA(iPred).multComp = multcompare(stats.ANOVA(iPred).anovaStats,'display','off');
end

%%
if plotFlag
cfg.clustCl      = [60 179 113; 145 203 55; 230 186 0; 150 110 70; 10 35 140; 90 115 210]./255;
figure;
subplot(1,5,[1 2]);
axs        = gca;
leafOrder  = optimalleaforder(stats.lnk,stats.distMat);
dh         = dendrogram(stats.lnk,'colorthreshold',0.35,'reorder',leafOrder,'labels',{}); % 'orientation','left',

% change colors
hold(axs, 'on')
lineColors = cell2mat(get(dh,'Color'));
colorList  = unique(lineColors, 'rows');

for color = 2:size(colorList,1)
  idx                = ismember(lineColors, colorList(color,:), 'rows');
  lineColors(idx, :) = repmat(cfg.clustCl(color-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for line = 1:numel(dh)
   dh(line).Color     = lineColors(line,:);
   dh(line).LineWidth = 1;
end
set(axs,'xtick',[],'xcolor','w')
ylabel('Eucl. dist.')


subplot(1,5,3)
xoffset  = 60;
ctBrain  = imread(analysisParams.ctBrainPath);
ctBrain  = repmat(ctBrain(:,xoffset+1:size(ctBrain,2)-xoffset,1),[1 1 3]);
ctBrain  = ctBrain+20;
ctBrain(ctBrain>255) = 255;
pxl(:,1) = (ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;
pxl(:,2) = -(AP.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(1);
pxl(:,3) = (-ML.*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2) - xoffset;

hold on
imshow(ctBrain)
for iK = 1:stats.nK
  plot(pxl(stats.clustID==iK,1),pxl(stats.clustID==iK,2),...
       '.','markersize',10,'color',cfg.clustCl(iK,:))
  plot(pxl(stats.clustID==iK,3),pxl(stats.clustID==iK,2),...
       '.','markersize',10,'color',cfg.clustCl(iK,:))
end
axis off; axis ij; axis image


subplot(1,5,[4 5]); hold on; axs=gca;
xt = zeros(1,size(stats.dataMat,2));
for iPt = 1:size(stats.dataMat,2)
  xt(iPt) = (iPt-1)*stats.nK + ceil(size(stats.dataMat,2)/2) + (iPt-1);
  for iK = 1:stats.nK
    x = (iPt-1)*stats.nK + iK + (iPt-1);
    if iPt == 1
      lh(iK)     = bar(x,100*mean(stats.dataMat(clustID==iK,iPt)),'facecolor',cfg.clustCl(iK,:),'edgecolor',cfg.clustCl(iK,:));
      lh_lbl{iK} = ['clust ' num2str(iK)];
    else
      bar(x,100*mean(stats.dataMat(clustID==iK,iPt)),'facecolor',cfg.clustCl(iK,:),'edgecolor',cfg.clustCl(iK,:))
    end
    errorbar(x,100*mean(stats.dataMat(clustID==iK,iPt)),std(100*stats.dataMat(clustID==iK,iPt))./sqrt(sum(clustID==iK)-1),...
             '-','color',cfg.clustCl(iK,:))
  end
  yl = get(axs,'ylim');
  text(xt(iPt),yl(2)*.98,sprintf('P = %1.2g',stats.ANOVA(iPred).pval)) % print p val
end

legend(lh,lh_lbl,'location','best'); legend('boxoff')
set(axs, 'xtick', xt, 'xticklabel', stats.labels)
xlim([0 xt(end)+ceil(stats.nK/2)+1])
ylabel(axs, 'Effect size (%)')
rotateXLabels(axs,30)
end
