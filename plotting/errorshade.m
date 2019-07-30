function errorshade(data,sem,linec,shadec,xaxis,lw)

% errorshade(data,sem,linec,shadec,xaxis,lw)
% plots data with error shades determined by sem. xaxis is optional, 
% default is vector going from 1 to length(data). colors should be a string
% following matlab conventions, or a 3-element vector (linec determines 
% line color, shadec determines shade color), defaults [0 0 0] and [.7 .7
% .7], respectively. lw is line width of th eplot, default 2
% LP, april 2011

% handle input arguments and defaults
if nargin < 2
    error('I need at least data and sem vectors')
elseif nargin == 2
    linec = 'k'; shadec = [.7 .7 .7]; xaxis = 1:length(data); lw = 2;
elseif nargin == 3
    if isstring(linec)
        error('please specify shade color')
    else shadec = [.7 .7 .7] + linec; shadec(shadec > 1) = 1;
        xaxis = 1:length(data);
    end
    lw = 2;
elseif nargin == 4
    xaxis = 1:length(data);
    lw = 2;
elseif nargin == 5
    lw = 2;
end

% make sure dimensions are right
if size(data,1) > size(data,2)
    data = data'; sem = sem';
end

% calculate shade polygons
fill_xvec = [xaxis fliplr(xaxis)];
fill_yvec = [data-sem fliplr(data+sem)]; 

% plot
fill(fill_xvec,fill_yvec,shadec,'lineStyle','none'); hold on
plot(xaxis,data,'Color',linec,'LineWidth',lw); %hold off
