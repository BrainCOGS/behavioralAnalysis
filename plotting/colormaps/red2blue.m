%RED to blue colormap

function map = red2blue(varargin)

rmap = [0:.01:1]'; omap = ones(size(rmap));
map  = [omap rmap rmap;flipud(rmap) flipud(rmap) omap];


map = colormap_helper(map, varargin{:});
