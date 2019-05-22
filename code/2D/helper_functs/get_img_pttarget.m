function [ f ] = get_img_pttarget(N, numTargets, width, changePct)
% N: dimension of image
% numTargets: number of targets to include
% changePct: rough percentage of targets to delete/move (0 is no changes)

f = zeros([N, N]);

gridLength = floor(N/ceil(sqrt(numTargets)));

x=-round(gridLength/2); y=-x;
for i=1:numTargets
    % increment to next point
    x = x + gridLength;
    if x > N
        x = round(gridLength/2);
        y = y + gridLength;
        if y > N
            break
        end
    end
    if changePct > rand()
        % make a change to whether the target is on the grid
        dx = randi([-5,5]); dy = randi([-5,5]);
        if 0.2 > rand()
            continue % don't put in target
        else
            f((x+dx):(x+dx+width-1), (y+dy):(y+dy+width-1)) = 1;
        end
    else
        f(x:(x+width-1),y:(y+width-1)) = 1;
    end
end

end

