function f1 = plotModelComparisonH(log_ML_array,log_OF_array,log_BFL_array,varargin)
    % Plot model-ranking metrics (ML, OF, BFL) in a grouped visual summary.
    % Normalizes scores by the best ML model to aid side-by-side comparison.
    %
    % Inputs:
    %   log_ML_array  - Log marginal likelihood values per model.
    %   log_OF_array  - Log objective-function scores per model.
    %   log_BFL_array - Log Bayes free-energy lower bound per model.
    %   varargin      - Optional plot controls (colors, labels, font, size).
    %
    % Outputs:
    %   f1 - Handle to the generated figure.

    % Set up:
    N_mod = length(log_ML_array);

    % Parse input arguments
    defaultColours = [colourLibrary(3,'gradient',[1 8]); colourLibrary(2)];
    defaultModelLabels = compose('Model %d',(1:N_mod));
    defaultFontSize = 14;
    defaultFontName = 'Times';
    defaultUnits = 'normalized';
    defaultWidth = 1;
    defaultHeight = 0.8;
    defaultFigNum = 3033;
    p = inputParser;
    addRequired(p,'log_ML_array');
    addRequired(p,'log_OF_array');
    addRequired(p,'log_BFL_array');
    addOptional(p,'colors',defaultColours);
    addOptional(p,'modelLabels',defaultModelLabels)
    addOptional(p,'fontSize',defaultFontSize);
    addOptional(p,'fontName',defaultFontName)
    addOptional(p,'figWidth',defaultWidth);
    addOptional(p,'figHeight',defaultHeight);
    addOptional(p,'figUnits',defaultUnits)
    addOptional(p,'figNum',defaultFigNum);

    parse(p,log_ML_array,log_OF_array,log_BFL_array,varargin{:});

    % Unpack optional arguments:
    colours = p.Results.colors;
    modelLabels = p.Results.modelLabels;
    fontSize = p.Results.fontSize;
    figW = p.Results.figWidth;
    figH = p.Results.figHeight;
    figU = p.Results.figUnits;
    figN = p.Results.figNum;

    % Make sure inputs are column vectors
    log_ML_array = log_ML_array(:);
    log_BFL_array = log_BFL_array(:);
    log_OF_array = log_OF_array(:);

    % Measure BFL and ML from the ML of the most likely model:
    [~,indML] = max(log_ML_array);
    log_BFL_array = log_BFL_array - log_ML_array(indML);
    log_ML_array = log_ML_array - log_ML_array(indML);

    % Prepare the labels for the plot
    labelArray = [log_BFL_array,log_OF_array,log_ML_array];

    f1 = figure(figN);
    set(f1,'units',figU,'position',[0 0 figW figH]);
    hold on
    
    % Set the bar width and spacing
    w = 0.28;
    buffer = 0.08;
    bar((1:N_mod)+w+buffer/2,log_ML_array,w,'FaceColor',colours(1,:));
    bar((1:N_mod)-w-buffer/2,log_BFL_array,w,'FaceColor',colours(3,:));
    bar(indML+w+buffer/2,log_ML_array(indML),w,'FaceColor',colours(1,:),'lineWidth',2);

    % Calculate the start and end points of each bar
    YEndPoints = [log_BFL_array,log_OF_array,log_ML_array];
    YEndPoints(:,2) = log_ML_array; % The OF bar should end at the ML tip
    YStartPoints = zeros(size(YEndPoints));
    YStartPoints(:,2) = log_BFL_array; % The OF bar should stem from the BFL tip
    % Correct the start and end points for the OF bar
    indBFLPositive = log_BFL_array >= 0;
    YEndPoints(indBFLPositive,2) = log_BFL_array(indBFLPositive);
    YEndPoints(~indBFLPositive,2) = log_ML_array(~indBFLPositive);
    YStartPoints(indBFLPositive,2) = log_ML_array(indBFLPositive);
    YStartPoints(~indBFLPositive,2) = log_BFL_array(~indBFLPositive);

    drawnow

    % Manually draw b(2) as rectangles
    for ii = 1:N_mod
        y_start = log_ML_array(ii); 
        y_end = log_BFL_array(ii); 
    
        % Draw rectangle (bar equivalent)
        rectangle('Position', [ii - w/2 , y_start, w, y_end - y_start], ...
                  'FaceColor', colours(2,:), 'EdgeColor', 'k'); 

        % Draw bounding lines
        plot([ii + 3*w/2 + buffer/2 , ii - w/2],[y_start y_start],'k')
        plot([ii + w/2 , ii - 3*w/2 - buffer/2],[y_end y_end],'k')
    end
                
    xlim([0.5 N_mod+0.5])

    ypadding = 5e-3*abs(diff(get(gca,'YLim')));
    
    hold on
    % Remove all boarders
    box off

    drawnow
    
    % Add labels to each of the bars
    for ii = 1:3
        ytips = YEndPoints(:,ii);
        xtips = (1:N_mod) + (ii-2).*(w + buffer/2);
        barHeight = YEndPoints(:,ii) - YStartPoints(:,ii);
        labels = compose('%.1f',labelArray(:,ii));
        t = text(xtips,ytips,labels,'VerticalAlignment','middle','FontSize',fontSize,'Rotation',90);
        textHeight = vertcat(t.Extent); textHeight = textHeight(:,4);

        %Check if each label fits within the bars. move it in if it
        %does, move it out if it doesn't.
        for jj = 1:length(t)
            pos = t(jj).Position;
            TextFitsInBar = abs(barHeight(jj)) > textHeight(jj) + ypadding;
            if TextFitsInBar
                if ytips(jj) < 0
                    ypos = ytips(jj) + ypadding;
                else
                    ypos = ytips(jj) - textHeight(jj) - ypadding;
                end
                pos(2) = ypos;
                set(t(jj),'Position',pos);
                set(t(jj),'Color','w')
            else
                if (ytips(jj) < 0)
                    ypos = ytips(jj) - textHeight(jj) - ypadding;
                else
                    ypos = ytips(jj) + ypadding;
                end
                pos(2) = ypos;
                set(t(jj),'Position',pos);
            end
        end

        % Add model label
        if ii == 2 % occam factor plot

            ylims = get(gca,'YLim');
            
            t2 = text(xtips,ylims(1)*ones(size(ytips)) + ypadding,modelLabels,'VerticalAlignment','middle','FontSize',fontSize,'FontWeight','bold');
            
            labelExtent = vertcat(t2.Extent);
            labelWidth = labelExtent(:,3);
            labelHeight = labelExtent(:,4);
            for jj = 1:length(t2)
                labelPos = t2(jj).Position;
                labelPos(1) = labelExtent(jj,1) - 0.5*labelWidth(jj);
                labelPos(2) = labelExtent(jj,2) - 0.5*labelHeight(jj);
                set(t2(jj),'Position',labelPos);
            end

            for jj = 1:length(t)
                pos = t(jj).Position;
                TextFitsInBar = abs(barHeight(jj)) > textHeight(jj) + ypadding;
                TextFitsInGap = abs(YStartPoints(jj,ii)) > textHeight(jj) + 2*ypadding;
                if (~TextFitsInBar) && (TextFitsInGap)
                    ypos = ytips(jj) - barHeight(jj) - ypadding;
                    pos(2) = ypos;
                    set(t(jj),'Position',pos);
                end
            end
        end

    end

    drawnow

    h = gca;
    set(h,'XColor','none')
    h.YAxis.Label.Color = [0,0,0];
    h.YAxis.TickLabelColor = [0 0 0];
    h.XAxis.Visible = 'off';
    yTicks = h.YAxis.TickValues;
    set(h,'FontSize',fontSize);
    % Add solid line at x=0
    yln = plot(xlim(h),[0,0],'k');
    % Add dashed line at other x ticks
    uistack(yln,'bottom');
    for ii = 1:length(yTicks)
        yln = plot(xlim(h),[yTicks(ii),yTicks(ii)],'--','color',[0.7 0.7 0.7]);
        uistack(yln,'bottom');
    end
    
    ax = gca;
    ax.Units = 'centimeters';

    t = ylabel('log(ML)','color',colours(1,:),'FontSize',fontSize);
    pos = t.Position;
    textWidth = vertcat(t.Extent); textWidth = 1.1*textWidth(3);
    text(pos(1)-2*textWidth,pos(2),'log(OF)','VerticalAlignment','top','HorizontalAlignment','center','color',colours(2,:),'FontSize',fontSize,'Rotation',90);
    text(pos(1)-3*textWidth,pos(2),'log(BFL)','VerticalAlignment','top','HorizontalAlignment','center','color',colours(3,:),'FontSize',fontSize,'Rotation',90);

    drawnow

    % Draw a line at the maximum marginal likelihood value:
    yln = plot(xlim(h),[max(log_ML_array),max(log_ML_array)],'color',colours(4,:));
    uistack(yln,'bottom');

    drawnow
end
