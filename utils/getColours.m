function rgb = getColours(colIdx,gradCnt)

    colours = [
             41  52 122; % 1. Blue
            253 129  83; % 2. Orange
            221  48  37; % 3. Red
            176 185 241; % 4. Blue2
            255 255 255; % 5. White
    ]/255;
    
    if nargin < 2
        gradCnt = 1;
    end
    
    if gradCnt == 1
        rgb = colours(colIdx,:);
    elseif gradCnt > 1 && ~isscalar(colIdx)
        M = length(colIdx);
        gradPoints = colours(colIdx,:);
        % Define the positions of the given colors along the gradient
        positions = linspace(1, gradCnt, M);
    
        % Interpolate each color channel separately
        rgb = zeros(gradCnt, 3);
        for channel = 1:3
            rgb(:, channel) = interp1(positions, gradPoints(:, channel), 1:gradCnt, 'linear');
        end
    else
        % Create a gradient from the seed colour to near white
        start = colours(colIdx,:);
        finish = [1,1,1];
        position = linspace(0,0.6,gradCnt).';
        rgb = start + position.*(finish-start);
    end

end