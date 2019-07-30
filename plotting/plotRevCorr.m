function plotRevCorr(lsrPerf,ca,lsrCl,printP)

% plotRevCorr(lsrPerf,ca,lsrCl,lsrSh,printP)
% plots laser ON and laser OFF logistic regression
%   lsrPerf is output structure of analyzeConcatLsrLog(). 
%           Must provide single cell entry in case of multi-area cell array
%   ca is axis handle
%   laserCl is color to plot laser data
%   laserSh is color to plot best-fitting curve
%   printP is boolean to print asterisks denoting p-value vs control

if nargin < 2; f  = figure; ca = gca(f);      end
if nargin < 3; lsrCl  = analysisParams.lsrCl; end
if nargin < 4; printP = true;                 end

try
  axes(ca); 
catch
  figure(ca);
end

hold on

xaxis1 = toBinCenters(lsrPerf.ctrl.logisticReg.bins);
plot([0 200],[0 0],'--','color',[.3 .3 .3],'linewidth',.25)
errorbar(xaxis1,lsrPerf.ctrl.logisticReg.values,lsrPerf.ctrl.logisticReg.sem,'-','color',analysisParams.ctrlCl,'linewidth',1);
errorbar(xaxis1,lsrPerf.lsr.logisticReg.values,lsrPerf.lsr.logisticReg.sem,'-','color',lsrCl,'linewidth',1);

if isfield(lsrPerf.stats,'logisticReg_pvals') && printP
  yl    = get(gca,'ylim');
  pvals = lsrPerf.stats.logisticReg_pvals';
  isSig = FDR(pvals);
  for iBin = 1:numel(pvals)
    if ~isSig(iBin)
      ptxt = '';%'n.s.';
    elseif isSig(iBin) && pvals(iBin) >= 0.01
      ptxt = '*';
    elseif  isSig(iBin) && pvals(iBin) >= 0.001 && pvals(iBin) < 0.01
      ptxt = '**';
    elseif isSig(iBin) && pvals(iBin) < 0.001
      ptxt = '***';
    end
    text(xaxis1(iBin),yl(2),ptxt,'fontsize',12,'horizontalAlignment','center')
  end
end

box off; set(gca,'fontsize',12,'xcolor','k','ycolor','k')
xlim([0 200]); %ylim([-.1 .35])
xlabel('y position (cm)','fontsize',14,'color','k')
ylabel('Weight on decision(a.u.)','fontsize',14,'color','k')