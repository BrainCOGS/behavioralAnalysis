function plotRevCorr_ctrl(rc,ca,cl)

% plotRevCorr_ctrl(rc,ca,cl)
% rc is the output structure of revCorr_logistic()
% ca is axis handle
% cl is color for data plotting

if nargin < 2 || isempty(ca)
    f  = figure;
    ca = gca(f);
end

if nargin < 3
  cl = analysisParams.ctrlCl;
end
try
    axes(ca); 
catch
    figure(ca);
end

hold on

xaxis1 = rc.bins(1:end-1) + mode(diff(rc.bins))/2;
errorbar(xaxis1,rc.values,rc.sem,'-','color',cl);

xlim([0 200]); 
box off; set(gca,'fontsize',12)
xlabel('Cue y (cm)','fontsize',14)
ylabel('Weight on decision','fontsize',14)