function [targPos,data] = gen_2D_pts(rand_targ,spacing,data)


fprintf('Generating targets...\n');
% Define the point targets here

data.maxR = data.maxWr/2 - spacing;              % Maximum range value to display (m)
data.minR = -data.maxWr/2 + spacing;             % Minimum range value to display (m)
data.maxX = data.maxWx/2 - spacing;              % Maximum x-range value to display (m)
data.minX = -data.maxWx/2 + spacing;             % Minimum x-range value to display (m)

targPos = zeros(data.N_targets,3);

if rand_targ
    x = data.minR + (data.maxR-data.minR).*rand(data.N_targets,1);
    y = data.minX + (data.maxX-data.minX).*rand(data.N_targets,1);
    
    vecii = x(1);
    vecjj = y(1);
    
    cnt = 2;
    for k = 2:data.N_targets
        thisX = x(k);
        thisY = y(k);
        distance = sqrt((thisX-vecii).^2 + (thisY-vecjj).^2);
        minDist = min(distance);
        
        iter = 0;
        while minDist < spacing
            
            thisX = data.minR + (data.maxR-data.minR).*rand;
            thisY = data.minX + (data.maxX-data.minX).*rand;
            
            distance = sqrt((thisX-vecii).^2 + (thisY-vecjj).^2);
            minDist = min(distance);
            
            iter = iter + 1;
            
            if iter > 500
                error('Decrease the spacing between targets');
            end
            
        end
        
        vecii(cnt) = thisX;
        vecjj(cnt) = thisY;
        cnt = cnt + 1;
        
    end
    
else
    vecii = [data.minR:spacing:data.maxR];%linspace(data.minR+10,data.maxR-10,num_rows)';
    vecjj = [data.minX:spacing:data.maxX]; %linspace(data.minX+10,data.maxX-10,num_rows)';
    
    [tmpx,tmpy] = meshgrid(vecii,vecjj);
    
    fprintf('Desired number of targets =  %d \n', data.N_targets);
    fprintf('Number of equally spaced targets that fit within the space = %d\n',...
        length(tmpx(:)));
    
    prompt = 'Change number of targets? [y/n]\n';
    increase = input(prompt,'s');
    
    if increase == 'y'
        
        targPos = zeros(length(tmpx(:)),3);
        
        vecii = tmpx(:);
        vecjj = tmpy(:);
        data.N_targets = length(tmpx(:));
        
    elseif increase == 'n'
        
        targPos = zeros(data.N_targets,3);
        tmp = tmpx(:);
        vecii = tmp(1:data.N_targets);
        tmp = tmpy(:);
        vecjj = tmp(1:data.N_targets);
        
    else
        sprintf('invalid answer \n');
        
    end
    
end

targPos(:,1) = vecii;
targPos(:,2) = vecjj;
targPos(:,3) = 0;

data.mag = randi(10,size(targPos,1),1);
data.amp = ones(size(targPos,1),1);
if data.iso == 0
    data.amp = data.amp.*data.mag; % for non isotropic points.
end

% Show the flight path and the targets
prompt = 'Display original target layout and flight path? [y/n]\n';
disp = input(prompt,'s');

if disp == 'y'
    
    figure
    subplot(1,3,1)
    plot(data.AntX,data.AntY,'r-');
    hold on
    for ii = 1:data.N_targets
        plot(targPos(ii,1),targPos(ii,2),'bx');
    end
    
    xlabel('x (m) - Range');
    ylabel('y (m) - Cross Range');
    
    X1 = [data.maxWr/2 data.maxWr/2 -data.maxWr/2 -data.maxWr/2 data.maxWr/2];
    X2 = [data.maxR data.maxR data.minR data.minR data.maxR];
    Y1 = [data.maxWx/2 -data.maxWx/2 -data.maxWx/2 data.maxWx/2 data.maxWx/2];
    Y2 = [data.maxX data.minX data.minX data.maxX data.maxX];
    
    plot(X1,Y1,'r'); % alias free scene
    plot(X2,Y2,'b'); % display window
    title('points and flight path')
    
    subplot(1,3,2);
    hold on;
    for ii = 1:data.N_targets
        plot(targPos(ii,1),targPos(ii,2),'bx');
    end
    plot(X1,Y1,'r');% alias free scene
    plot(X2,Y2,'b'); % display window
    
    title('points only')
    
    subplot(1,3,3);
    plot(data.AntX,data.AntY,'r-o');
    title('flight path')
    
    % aliasing check
    if (sum(abs(X2) >= abs(X1))> 0 || sum(abs(Y2) >= abs(Y1))>0)
        error('Error: aliasing. Consider decreasing N_targets or spacing.');
    end
    
    drawnow
end

end

