function [targPos,data] = gen_3D_pts(rand_targ,spacing,data)


fprintf('Generating targets...\n');
% Define the point targets here
if mod(sqrt(data.N_targets),2) == 0
    tmp = sqrt(data.N_targets)/2;
    left_loc = -spacing*(tmp-1);
    right_loc = spacing*tmp;
else
    tmp = floor(sqrt(data.N_targets)/2);
    left_loc = -spacing*tmp;
    right_loc = spacing*tmp;
end

data.maxR = data.maxWr/2 - spacing;              % Maximum range value to display (m)
data.minR = -data.maxWr/2 + spacing;             % Minimum range value to display (m)
data.maxX = data.maxWx/2 - spacing;              % Maximum x-range value to display (m)
data.minX = -data.maxWx/2 + spacing;             % Minimum x-range value to display (m)
data.maxZ = right_loc + spacing;
data.minZ = left_loc - spacing;

targPos = zeros(data.N_targets,3);

if rand_targ
    
    fprintf('Spacing = %d\n',spacing); 
    
    x = data.minR + (data.maxR-data.minR).*rand(data.N_targets,1);
    y = data.minX + (data.maxX-data.minX).*rand(data.N_targets,1);
    z = data.minX + (data.maxZ-data.minZ).*rand(data.N_targets,1);
    
    vecii = x(1);
    vecjj = y(1);
    veckk = z(1);
    
    cnt = 2;
    for k = 2:data.N_targets
        
        thisX = x(k);
        thisY = y(k);
        thisZ = z(k);
        
        distance = sqrt((thisX-vecii).^2 + (thisY-vecjj).^2 + ...
            (thisZ - veckk).^2);
        
        minDist = min(distance);
        
        iter = 0;
        while minDist < spacing
            
            thisX = data.minR + (data.maxR-data.minR).*rand;
            thisY = data.minX + (data.maxX-data.minX).*rand;
            thisZ = data.minZ + (data.maxZ-data.minZ).*rand;
            
            distance = sqrt((thisX-vecii).^2 + (thisY-vecjj).^2 + ...
                (thisZ - veckk).^2);
            
            minDist = min(distance);
            
            iter = iter + 1;
            
            if iter > 500
                error('Decrease your spacing between targets');
            end
            
        end
        
        vecii(cnt) = thisX;
        vecjj(cnt) = thisY;
        veckk(cnt) = thisZ;
        
        fprintf('The distance between point %d and all others = %.4f\n',...
            k,minDist);
        
        cnt = cnt + 1;
        
    end
    
else
    vecii = [data.minR:spacing:data.maxR];%linspace(data.minR+10,data.maxR-10,num_rows)';
    vecjj = [data.minX:spacing:data.maxX]; %linspace(data.minX+10,data.maxX-10,num_rows)';
    veckk = [data.minZ:spacing:data.maxZ];
    
    [tmpx,tmpy,tmpz] = meshgrid(vecii,vecjj,veckk);
    
    fprintf('Desired number of targets =  %d \n', data.N_targets);
    fprintf('Number of equally spaced targets that fit within the space = %d\n',...
        length(tmpx(:)));
    
    prompt = 'increase number of targets? [y/n]\n';
    increase = input(prompt,'s');
    
    if increase == 'y'
        
        targPos = zeros(length(tmpx(:)),3);
        
        vecii = tmpx(:);
        vecjj = tmpy(:);
        veckk = tmpz(:);
        data.N_targets = length(tmpx(:));
        
    elseif increase == 'n'
        
        targPos = zeros(data.N_targets,3);
        tmp = tmpx(:);
        vecii = tmp(1:data.N_targets);
        tmp = tmpy(:);
        vecjj = tmp(1:data.N_targets);
        tmp = tmpz(:);
        veckk = tmp(1:data.N_targets);
        
    else
        sprintf('invalid answer \n');
        
    end
end

targPos(:,1) =vecii;
targPos(:,2) = vecjj;
targPos(:,3) = veckk;

data.mag = randi(10,size(targPos,1),1);
data.amp = ones(size(targPos,1),1);
if data.iso == 0
    data.amp = data.amp.*data.mag; % for non isotropic points.
end

X1 = [data.maxWr/2 data.maxWr/2 -data.maxWr/2 -data.maxWr/2 data.maxWr/2];
X2 = [data.maxR data.maxR data.minR data.minR data.maxR];
Y1 = [data.maxWx/2 -data.maxWx/2 -data.maxWx/2 data.maxWx/2 data.maxWx/2];
Y2 = [data.maxX data.minX data.minX data.maxX data.maxX];

% aliasing check
if (sum(abs(X2) >= abs(X1))> 0 || sum(abs(Y2) >= abs(Y1))>0)
    error('Error: aliasing in xy-plane. Consider decreasing N_targets or spacing.');
end

% Show the flight path and the targets
prompt = 'Display original target layout and flight path? [y/n]\n';
disp = input(prompt,'s');

if disp == 'y'
    figure;
    
    subplot(1,3,1);
    plot3(data.AntX,data.AntY,data.AntZ,'r-'); hold on;
    scatter3(targPos(:,1),targPos(:,2),...
        targPos(:,3),100,data.amp,'.')
    grid('on')
    title('points and flight path')
    
    subplot(1,3,2);
    scatter3(targPos(:,1),targPos(:,2),...
        targPos(:,3),100,data.amp,'.')
    title('points only')
    
    subplot(1,3,3);
    plot3(data.AntX,data.AntY,data.AntZ,'r-o');
    grid('on')
    title('flight path')
    drawnow;
    
end

end

