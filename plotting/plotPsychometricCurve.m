function plotPsychometricCurve(lsrPerf,ca,lsrCl,lsrSh,printP)

% plotPsychometricCurve(lsrPerf,ca,lsrCl,lsrSh,printP)
% plots laser ON and laser OFF psychometrics
%   lsrPerf is output structure of analyzeConcatLsrLog(). 
%           Must provide single cell entry in case of multi-area cell array
%   ca is axis handle
%   laserCl is color to plot laser data
%   laserSh is color to plot best-fitting curve
%   printP is boolean to print text w/ p-value comparing the slopes

%%
if nargin < 2; figure; ca = gca;                end
if nargin < 3; lsrCl = analysisParams.lsrCl;    end
if nargin < 4; lsrSh = analysisParams.lsrShade; end
if nargin < 5; printP = true;                   end

%%
try
  axes(ca);
catch
  figure(ca);
end

hold on

plot([-15 15],[50 50],':','color',[.7 .7 .7])
plot([0 0],[0 100],':','color',[.7 .7 .7])

x = lsrPerf.ctrl.psychometric.perfPsych_xaxisBins;
y = 100.*lsrPerf.ctrl.psychometric.perfPsychJ_binned;
l = 100.*-lsrPerf.ctrl.psychometric.perfPsychJ_binnedSEM(:,1)+y; 
u = 100.*lsrPerf.ctrl.psychometric.perfPsychJ_binnedSEM(:,2)-y; 

errorbar(x,y,l,u,'o','color',analysisParams.ctrlCl,'markersize',6,'markerfacecolor',analysisParams.ctrlCl)
plot(lsrPerf.ctrl.psychometric.fit.xaxis,100.*lsrPerf.ctrl.psychometric.fit.curve,'-','color',analysisParams.ctrlShade,'linewidth',1)

x = lsrPerf.lsr.psychometric.perfPsych_xaxisBins;
y = 100.*lsrPerf.lsr.psychometric.perfPsychJ_binned;
l = 100.*-lsrPerf.lsr.psychometric.perfPsychJ_binnedSEM(:,1)+y; 
u = 100.*lsrPerf.lsr.psychometric.perfPsychJ_binnedSEM(:,2)-y; 

errorbar(x,y,l,u,'o','color',lsrCl,'markersize',6,'markerfacecolor',lsrCl)
plot(lsrPerf.lsr.psychometric.fit.xaxis,100.*lsrPerf.lsr.psychometric.fit.curve,'-','color',lsrSh,'linewidth',1)

box off; set(gca,'fontsize',12,'xcolor','k','ycolor','k','ytick',0:25:100)
xlim([-lsrPerf.ctrl.psychometric.perfPsych_xaxisBins(end)-1 ...
    lsrPerf.ctrl.psychometric.perfPsych_xaxisBins(end)+1]);
ylim([0 100])
xlabel('\Delta towers','fontsize',14,'color','k')
ylabel('Went right (%)','fontsize',14,'color','k')

text(-13,90,'Laser off','color',analysisParams.ctrlCl,'fontsize',12)
text(-13,83,'Laser on','color',lsrCl,'fontsize',12)

if isfield(lsrPerf.stats.bootPerf,'p_slope') && printP
  text(13,20,sprintf('p(slope) = %1.2g',lsrPerf.stats.bootPerf.p_slope),...
       'fontsize',11,'horizontalAlignment','right')
end

if isfield(lsrPerf,'cueHalf')
  x = lsrPerf.cueHalf.psychometric.perfPsych_xaxisBins;
  y = 100.*lsrPerf.cueHalf.psychometric.perfPsychJ_binned;
  l = 100.*-lsrPerf.cueHalf.psychometric.perfPsychJ_binnedSEM(:,1)+y; 
  u = 100.*lsrPerf.cueHalf.psychometric.perfPsychJ_binnedSEM(:,2)-y; 

  errorbar(x,y,l,u,'o','color',analysisParams.cueHalfCl,'markersize',6,...
                   'markerfacecolor',analysisParams.cueHalfCl)
  plot(lsrPerf.cueHalf.psychometric.fit.xaxis,100.*lsrPerf.cueHalf.psychometric.fit.curve,...
                   '-','color',analysisParams.cueHalfShade,'linewidth',1)
  text(-13,76,'Deleted towers','color',analysisParams.cueHalfCl,'fontsize',12)
end
