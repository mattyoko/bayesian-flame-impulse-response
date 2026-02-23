function cornerHeatmap(X, bMAP, CbMAP, paramNames, nbins)
% Plot pairwise posterior marginals as a corner heatmap.
%
% Inputs:
%   X          - Sample matrix of size [N x D].
%   bMAP       - MAP parameter vector of length D.
%   CbMAP      - MAP covariance matrix of size [D x D].
%   paramNames - Optional 1xD cell array of axis labels.
%   nbins      - Optional number of bins (default 40).

cmap = getColours([3,5],50);
cmap = flipud(cmap);

if nargin < 2 
    paramNames = compose('p_%d', 1:size(X,2));
end
if nargin < 3 
    nbins = 40; 
end

D = size(X,2);
tiledlayout(D, D, 'TileSpacing','none','Padding','compact');

for ii = 1:D
    for jj = 1:D
        ax = nexttile((ii-1)*D + jj);
        if ii == jj

            mu = bMAP(jj);
            sig = sqrt(CbMAP(jj,ii));

            edges = linspace(mu-4*sig,mu+4*sig,nbins);

            histogram(X(:,jj), edges, 'Normalization','pdf', ...
                'EdgeColor', 'none', ...
                'FaceColor', getColours(3),...
                'FaceAlpha', 0.5);
           
            hold on

            histogram(X(:,jj), edges, 'Normalization','pdf', ...
                'DisplayStyle','stairs', ...
                'EdgeColor', getColours(3), ...
                'EdgeAlpha', 1.0,...
                'LineWidth', 1.0);
            set(ax,'YTick',[]);
            
            b = linspace(mu-4*sig,mu+4*sig);
            pb = normpdf(b, mu, sig);
            plot(b, pb, 'color',getColours(1), 'LineWidth', 1.5);
            xlim([mu - 4*sig,mu + 4*sig])
            ylim([0 1.2*max(pb)])
            box on
        else
            mux  = bMAP(jj);
            muy  = bMAP(ii);
            sigx = sqrt(CbMAP(jj,jj));
            sigy = sqrt(CbMAP(ii,ii));

            xedges = linspace(mux - 4*sigx, mux + 4*sigx, nbins+1);
            yedges = linspace(muy - 4*sigy, muy + 4*sigy, nbins+1);

            counts = histcounts2(X(:,jj), X(:,ii), xedges, yedges, 'Normalization','pdf')';
            counts(~isfinite(counts)) = 0;

            % --- alpha map ---
            cmax = prctile(counts(:), 99);
            if ~isfinite(cmax) || cmax <= 0, cmax = max(counts(:)); end
            cmax = max(cmax, eps);

            t = counts ./ cmax;                
            t = max(0, min(1, t));

            alphaPower = 0.6;                  
            alphaMin   = 0.00;                 
            alphaMax   = 0.95;                  
            Aalpha = alphaMin + (alphaMax-alphaMin) * (t .^ alphaPower);
            Aalpha(t < 0.02) = 0;

            % --- draw as image ---
            hIm = imagesc(ax, [xedges(1) xedges(end)], [yedges(1) yedges(end)], counts);
            set(ax, 'YDir','normal');
            set(hIm, 'AlphaData', Aalpha, 'AlphaDataMapping','none');

            hold on

            % --- overlay Laplace ellipses ---
            mu = [mux; muy];
            C = CbMAP([jj,ii],[jj,ii]);
            drawCovEllipse(mu,C,1,getColours(1));
            drawCovEllipse(mu,C,2,getColours(1));
            drawCovEllipse(mu,C,3,getColours(1));

            xlim([xedges(1) xedges(end)]);
            ylim([yedges(1) yedges(end)]);
            set(ax,'color',cmap(1,:))
            box on
            grid off
        end
        if ii < D 
            set(ax,'XTickLabel',[]); 
        else 
            xlabel(paramNames{jj},'interpreter','latex'); 
            set(ax,'XTick',[mux - 2*sigx, mux, mux + 2*sigx])
            xtickformat('%.2f')
            xtickangle(45)
        end
        if jj > 1 
            set(ax,'YTickLabel',[]); 
        else 
            ylabel(paramNames{ii},'interpreter','latex'); 
            if ii > 1
                set(ax,'YTick',[muy - 2*sigy, muy, muy + 2*sigy])
                ytickformat('%.2f')
            end
        end
    end
end
colormap(cmap);
end

function drawCovEllipse(mu, C, k, col)
% Draw a k-sigma covariance ellipse for a 2D Gaussian.
%
% Inputs:
%   mu  - Mean vector [2 x 1].
%   C   - Covariance matrix [2 x 2].
%   k   - Sigma level for the ellipse radius.
%   col - Line color specification.
%
% Outputs:
%   None. The ellipse is plotted in the current axes.
    if any(~isfinite(mu)) || any(~isfinite(C(:)))
        return; 
    end
    [V, D] = eig((C+C')/2);
    D = max(D,0); 
    t = linspace(0, 2*pi, 200);
    a = k*sqrt(max(D(1,1),0)); b = k*sqrt(max(D(2,2),0));
    ellipse = V * [a*cos(t); b*sin(t)];
    x = mu(1) + ellipse(1,:); y = mu(2) + ellipse(2,:);
    plot(x, y, 'color', col, 'LineWidth', 1);
end