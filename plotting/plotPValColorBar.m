function c = plotPValColorBar(ca)

% c = plotPValColorBar(ca)
% funnction to plot stand-alone pvalue color bar for inactivation maps
% ca is axis handle, c is colorbar handle

if nargin < 1
    f        = figure;
    ca       = gca(f);
end


cmap    = flipud(colormap(analysisParams.colormap));

axes(ca)
imagesc([],[],[],log10(analysisParams.cbarLim))
colormap(cmap)
c = colorbar('Limits',log10(analysisParams.cbarLim),'location','SouthOutside',...
    'fontsize',11,'Position',[0.72 0.28 0.18 0.04]);
c.Label.String = 'log10(P-value)';
c.Label.FontSize = 12;
set(ca,'visible','off')
