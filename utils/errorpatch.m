function hl =  errorpatch(x,mean,lo,hi,varargin)
% Plot a mean line with asymmetric shaded uncertainty bounds.
% Additional line/patch styling can be passed through name-value pairs.
%
% Inputs:
%   x       - Horizontal coordinates.
%   mean    - Mean line values.
%   lo      - Lower deviation from mean.
%   hi      - Upper deviation from mean.
%   varargin- Optional style overrides:
%             * 'color' : base color
%             * 'Line'  : cell array of line property pairs
%             * 'Patch' : cell array of patch property pairs
%
% Outputs:
%   hl - Handle to the plotted mean line.

% Make sure all inputs are column vectors
x = x(:);
mean = mean(:);
lo = lo(:);
hi = hi(:);

% Create the line and patch elements
hp = patch([x ; flipud(x)] , [mean + hi ; flipud(mean - lo)],[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.2);
hold on
hl = plot(x,mean,'lineWidth',2);

% Set up input parser to manage 'Line' and 'Patch' property groups
p = inputParser;
addParameter(p, 'color' , 'b');
addParameter(p, 'Line', {}, @(x) iscell(x) && mod(length(x), 2) == 0);
addParameter(p, 'Patch', {}, @(x) iscell(x) && mod(length(x), 2) == 0);
parse(p, varargin{:});

% Extract parsed properties for line and patch
color = p.Results.color;
lineProps = p.Results.Line;
patchProps = p.Results.Patch;

% Separate properties into names and values for the line
linePropNames = lineProps(1:2:end);
linePropValues = lineProps(2:2:end);

% Separate properties into names and values for the patch
patchPropNames = patchProps(1:2:end);
patchPropValues = patchProps(2:2:end);


if ~isempty(color)
    set(hl,'color',color);
    set(hp,'FaceColor',color);
end
if ~isempty(linePropNames)
    set(hl,linePropNames,linePropValues);
end
if ~isempty(patchPropNames)
    set(hp,patchPropNames,patchPropValues);
end

end