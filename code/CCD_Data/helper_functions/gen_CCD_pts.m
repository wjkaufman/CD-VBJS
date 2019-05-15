function [targPos_CCD,data] = gen_CCD_pts(targPos,dimension,data)

targPos_CCD = targPos;
data.amp_CCD = data.amp;

% targets to be removed
ind = randperm(data.N_targets,data.remove);
targPos_CCD(ind,:) = [];
data.amp_CCD(ind) = [];

% targets to change location
[data.N_targets_CCD,~] = size(targPos_CCD);
ind = randperm(data.N_targets_CCD,data.num_ccd-data.remove);

move_x = min(targPos_CCD(:,1)) + (max(targPos_CCD(:,1))-...
    min(targPos_CCD(:,1)))*rand(data.num_ccd-data.remove,1);
move_y = min(targPos_CCD(:,2)) + (max(targPos_CCD(:,2))-...
    min(targPos_CCD(:,2)))*rand(data.num_ccd-data.remove,1);
if dimension == 3
    move_z = data.maxZ + (data.maxZ-...
        data.minZ)*rand(data.num_ccd-data.remove,1);
    targPos_CCD(ind,3) = targPos_CCD(ind,3) + move_z/2;
    ind_bad_z = find(targPos_CCD(:,3)>data.maxZ|targPos_CCD(:,3)<data.minZ);
else
    ind_bad_z = [];
end

targPos_CCD(ind,1) = targPos_CCD(ind,1) + move_x/2;
targPos_CCD(ind,2) = targPos_CCD(ind,2) + move_y/2;

% make sure targets still fall within the imaging grid
ind_bad_x = find(targPos_CCD(:,1)>data.maxR|targPos_CCD(:,1)<data.minR);
ind_bad_y = find(targPos_CCD(:,2)>data.maxX|targPos_CCD(:,2)<data.minX);

while ~isempty(ind_bad_x)
    targPos_CCD(ind_bad_x,1) = targPos_CCD(ind_bad_x,1)/2;
    ind_bad_x = find(targPos_CCD(:,1)>data.maxR|targPos_CCD(:,1)<data.minR);
end

while ~isempty(ind_bad_y)
    targPos_CCD(ind_bad_y,2) = targPos_CCD(ind_bad_y,2)/2;
    ind_bad_y = find(targPos_CCD(:,2)>data.maxX|targPos_CCD(:,2)<data.minX);
end

while ~isempty(ind_bad_z)
    targPos_CCD(ind_bad_z,3) = targPos_CCD(ind_bad_z,3)/2;
    ind_bad_z = find(targPos_CCD(:,3)>data.maxZ|targPos_CCD(:,3)<data.minZ);
end

% display taret locations
if dimension == 2
    
    X1 = [data.maxWr/2 data.maxWr/2 -data.maxWr/2 -data.maxWr/2 data.maxWr/2];
    X2 = [data.maxR data.maxR data.minR data.minR data.maxR];
    Y1 = [data.maxWx/2 -data.maxWx/2 -data.maxWx/2 data.maxWx/2 data.maxWx/2];
    Y2 = [data.maxX data.minX data.minX data.maxX data.maxX];
    
    figure;
    subplot(1,3,1);
    hold on;
    for ii = 1:data.N_targets
        plot(targPos(ii,1),targPos(ii,2),'bx');
    end
    plot(X1,Y1,'r');% alias free scene
    plot(X2,Y2,'b'); % display window
    
    title('Reference Scene')
    
    subplot(1,3,2);
    hold on;
    for ii = 1:data.N_targets_CCD
        plot(targPos_CCD(ii,1),targPos_CCD(ii,2),'gs');
    end
    plot(X1,Y1,'r');% alias free scene
    plot(X2,Y2,'b'); % display window
    title('Changed Scene')
    
    subplot(1,3,3);
    hold on;
    for ii = 1:data.N_targets
        plot(targPos(ii,1),targPos(ii,2),'bx');
    end
    
    for ii = 1:data.N_targets_CCD
        plot(targPos_CCD(ii,1),targPos_CCD(ii,2),'gs');
    end
    plot(X1,Y1,'r');% alias free scene
    plot(X2,Y2,'b'); % display window
    title('Comparison')
    
elseif dimension == 3
    
    figure;
    subplot(1,2,1);
    scatter3(targPos(:,1),targPos(:,2),...
        targPos(:,3),100,data.amp,'.')
    title('Reference Scene')
    
    subplot(1,2,2);
    scatter3(targPos_CCD(:,1),targPos_CCD(:,2),...
        targPos_CCD(:,3),100,data.amp_CCD,'.')
    title('Changed Scene')
    
end

end