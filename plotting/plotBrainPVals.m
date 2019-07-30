function [pmat,cbarHandle] = plotBrainPVals(ML,AP,pval,effectSize,alpha_correct,thFlag,type,ca,plotBar)

% [pmat,cbarHandle] = plotBrainPVals(ML,AP,pval,effectSize,alpha_correct,thFlag,type,ca,plotBar)
% Input:
%   ML list of ML coordinates in mm
%   AP list of AP coordinates in mm
%   pval list of p-values, same numel as ML and AP
%   effect size list of magnitude and sign of bias / % correct
%   alpha_correct is the alpha obtained by the FDR correction procedure
%   thFlag is flag to threshold colors at alpha_correct
%   type is 'es_size' for circle size proportional to effect size and color
%        to p-value; vice-versa for 'p_size'
%   ca is axis handle if desired
%   plotBar is to plot color bar
%
% Output
%   pmat RGB matrix with figure
%   cbarHandle is handle for color bar if applicable

%% defaults
if nargin < 4
  effectSize    = ones(size(pval));
end
if nargin < 5
  alpha_correct = analysisParams.bootalpha;
end
if nargin < 6
  thFlag        = 1;
end
if nargin < 7
  type          = 'es_size';% circ diameter denotes p or effectsize?
end
if nargin < 8
  f             = figure;
  ca            = gca(f);
end
if nargin < 9
  plotBar = true;
end

%% generate image
effectSize(effectSize > 1)  = 1;
effectSize(effectSize < -1) = -1;

xoffset = 60;
ctBrain = imread(analysisParams.ctBrainPath);
pmat    = repmat(ctBrain(:,xoffset+1:size(ctBrain,2)-xoffset,1),[1 1 3]);
pmat    = pmat+20;
pmat(pmat>255) = 255;

% colormap; invert so brightest is lowest P-val
if ~strcmpi(analysisParams.colormap,'red2blue')
  cmap    = flipud(colormap(analysisParams.colormap));
else
  cmap    = colormap(analysisParams.colormap);
end

cmap = cmap.*255;

switch type
  case 'es_size'
    if analysisParams.cbarLim(1) < 0
      pLookUp = linspace(-log10(abs(analysisParams.cbarLim(1))),log10(analysisParams.cbarLim(2)),length(cmap)); % p will be in log scale
    else
      pLookUp = linspace(log10(analysisParams.cbarLim(1)),log10(analysisParams.cbarLim(2)),length(cmap)); % p will be in log scale
    end
    
  case 'p_size'
    if isempty(effectSize)
      pLookUp = zeros(1,length(cmap));
    else
      if analysisParams.cbarLim(1) < 0
        pLookUp = linspace(-max(abs(effectSize)),max(abs(effectSize)),length(cmap)); % p will be in log scale
      else
        pLookUp = linspace(analysisParams.cbarLimEs(1),analysisParams.cbarLimEs(2),length(cmap)); % p will be in log scale
      end
    end
    
end

%% add color-coded p-value circles
for ii = 1:length(pval)
  centerPxl = getPxlCoord(ML(ii),AP(ii))-[xoffset 0];
  thisp     = pval(ii);
  if thFlag && strcmpi(type,'es_size')
    if thisp > alpha_correct
      thisp = 1;
    end
  end
  
  pmat      = drawPcolors(thisp,effectSize(ii),pmat,centerPxl,cmap,pLookUp,type);
end

%% reset colormap for colorbar
if ~strcmpi(analysisParams.colormap,'red2blue')
  cmap    = flipud(colormap(analysisParams.colormap));
else
  cmap    = colormap(analysisParams.colormap);
end

%% add p-value circle size caption
switch type
  case 'es_size'
    esCaption = [0.05 .25];
    MLcap     = [-6.2 -6.2]; %[-4.5 -3.5 -2.5];
    APcap     = [4.5 3.75];
    for ii = 1:length(esCaption)
      centerPxl = getPxlCoord(MLcap(ii),APcap(ii));
      pmat      = drawPcolors(0,esCaption(ii),pmat,centerPxl,zeros(size(cmap)),pLookUp,type);
    end
    
    % plot
    
    barLookUp = linspace(0,1,length(cmap));
    ticks = [0 barLookUp(find(pLookUp<=abs(log10(analysisParams.bootalpha)),1,'first')) ...
      .5 barLookUp(find(pLookUp>=log10(analysisParams.bootalpha),1,'last'))  1];
    tickLbls = {sprintf('1e%d',pLookUp(end)),sprintf('%1.3f',analysisParams.bootalpha),'1',sprintf('%1.3f',analysisParams.bootalpha),num2str(10^pLookUp(end))};
    
    axes(ca)
    colormap(cmap);
    imshow(imresize(pmat,2))
    
    if plotBar
      cbarHandle = smallcolorbar(gca,'southoutside');
      set(cbarHandle,'ticks',ticks,'tickLabels',tickLbls)
      cbarHandle.Label.String = '\Delta < 0       P-value       \Delta > 0'; %(- <- -> +) or (L <- -> R)
    else
      cbarHandle = [];
    end
    
    text(1,-45,'\Delta:','fontsize',12,'horizontalAlignment','left')
    text(80,20,'5  %','fontsize',12,'horizontalAlignment','left')
    text(80,100,'25 %','fontsize',12,'horizontalAlignment','left')
  case 'p_size'

    % plot
    barLookUp = linspace(0,1,length(cmap));
    ticks = [0 .25 .5 .75 1];
    tickLbls = {num2str(round(pLookUp(1)*100)),num2str(round(pLookUp(1)/2*100)),...
      '0',num2str(round(pLookUp(end)/2*100)),num2str(round(pLookUp(end)*100))};
    
    axes(ca)
    colormap(cmap);
    imshow(imresize(pmat,2))
    
    if plotBar
      cbarHandle = smallcolorbar(gca,'southoutside');
      set(cbarHandle,'ticks',ticks,'tickLabels',tickLbls)
      cbarHandle.Label.String = '\Delta (%)'; %(- <- -> +) or (L <- -> R)
      text(80,15,'\Delta:','fontsize',11,'horizontalAlignment','left')
      text(80,50,'5%','fontsize',11,'horizontalAlignment','left')
      text(80,100,'50%','fontsize',11,'horizontalAlignment','left')
    end
end


end

%% transform [ML AP] into pixel space
function pxl = getPxlCoord(ML,AP)

pxl(1) = round(ML*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(2);
pxl(2) = -round(AP*analysisParams.ctBrainPxlPerMM) + analysisParams.ctBrainBregma(1);

end

%% draw coordinate grid
function pmat = drawgrid(pmat)

AP = -4:4;
ML = -3:4;

for ii = 1:length(AP)
  pxl = getPxlCoord(0,AP(ii));
  pmat(:,pxl(2)) = .6;
end

for ii = 1:length(ML)
  pxl = getPxlCoord(ML(ii),0);
  pmat(pxl(1),:) = .6;
end

end

%% look up color map and plot color-coded p-value circle
function pmat = drawPcolors(p,es,pmat,centerPxl,cmap,pLookUp,type)
try
  switch type
    case 'es_size'
      if p < 10^pLookUp(end)
        p = 10^pLookUp(end);
      end
      
      plog = log10(p).*sign(es);
      if plog > 0
        cl   = cmap(find(pLookUp<=plog,1,'first'),:);
      elseif plog < 0
        cl   = cmap(find(pLookUp>=plog,1,'last'),:);
      elseif plog == 0
        cl   = [255 255 255];
      end
      rad  = round(analysisParams.radiusMM(2)*analysisParams.ctBrainPxlPerMM); % max radius
      rad  = round(rad.*abs(es));
      if rad<analysisParams.radiusMM(1)*analysisParams.ctBrainPxlPerMM
        rad = round(analysisParams.radiusMM(1)*analysisParams.ctBrainPxlPerMM); % min radius
      end
      
      for ii = centerPxl(2)-rad:centerPxl(2)+rad
        for jj = centerPxl(1)-rad:centerPxl(1)+rad
          if ((ii - centerPxl(2))^2 + (jj- centerPxl(1))^2) <= rad^2
            try pmat(max([1 ii-2]),max([1 jj+3]),:) = reshape(cl,[1 1 3]); catch; keyboard; end
          end
        end
      end
      
      
    case 'p_size'
      
      temp = linspace(-log10(abs(analysisParams.cbarLim(1))),log10(analysisParams.cbarLim(2)),length(cmap));
      if p < 10^temp(end)
        p = 10^temp(end);
      end
      
      plog = abs(log10(p));
      rad  = round(analysisParams.radiusMM(2)*analysisParams.ctBrainPxlPerMM); % max radius
      rad  = round(rad.*plog/10);
      if rad<analysisParams.radiusMM(1)*analysisParams.ctBrainPxlPerMM
        rad = round(analysisParams.radiusMM(1)*analysisParams.ctBrainPxlPerMM); % min radius
      end
      if rad > 23; keyboard; end
      cl   = cmap(find(pLookUp>=es,1,'first'),:);
      
      for ii = centerPxl(2)-rad:centerPxl(2)+rad
        for jj = centerPxl(1)-rad:centerPxl(1)+rad
          if ((ii - centerPxl(2))^2 + (jj- centerPxl(1))^2) <= rad^2
            pmat(max([1 ii-2]),max([1 jj+3]),:) = reshape(cl,[1 1 3]);
          end
        end
      end
  end
catch
  keyboard
end
end
