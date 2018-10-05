function [ f ] = get_img( funct,N)

if strcmp(funct,'phantom')
    orig_img = phantom(N);
elseif strcmp(funct,'multishape')
    orig_img = zeros(200,200);
    orig_img(10:160,20:25) = 1;
    orig_img(60:180,40:60) = .5;
    orig_img(120:140,80:180) = 1;
    orig_img(160:180,120:190) = .75;
    x = linspace(0,200,200);
    y = linspace(0,200,200);
    for i = 1:200
        for j = 1:200
            if ((x(i)-60)^2 +(y(j)-130)^2 < 40^2)
                orig_img(i,j) = 1;
            end
        end
    end
    orig_img = orig_img + 1;
elseif strcmp(funct,'circ_sq')
    % circle/square (on 10 to 20)
    x = linspace(-1,1,N);
    y = linspace(-1,1,N);
    
    orig_img = 10*ones(N,N);
    r = .75;
    
    for i = 1:N
        for j = 1:N
            
            if ((x(i))^2 + (y(j))^2 < r^2)
                orig_img(i,j) = orig_img(i,j) + 10;
            end
            
            if (x(i)<.25 & x(i) > -.25)&(y(j)< .25 & y(j) > -.25)
                orig_img(i,j) = orig_img(i,j) - 5;
            end
        end
    end
elseif strcmp(funct,'multicirc')
    x = linspace(-50,50,N);
    y = linspace(-50,50,N);
    
    orig_img = 10*ones(N,N);
    
    r1 = 10;
    r2 = 5;
    r3 = 20;
    r4 = 15;
    for i = 1:N
        for j = 1:N
            if ((x(i)+30)^2 + (y(j)-30)^2 < r1^2)
                orig_img(i,j) = orig_img(i,j) + 5;
            end
            
            if ((x(i)-20)^2 + (y(j)-20)^2 < r4^2)
                orig_img(i,j) = orig_img(i,j) + 10;
            end
            
            
            if ((x(i)+10)^2 + (y(j)+10)^2 < r3^2)
                orig_img(i,j) = orig_img(i,j) + 15;
            end
            
            if ((x(i)-30)^2 + (y(j)+30)^2 < r2^2)
                orig_img(i,j) = orig_img(i,j) + 20;
            end
            
        end
    end
    
elseif strcmp(funct,'hill_sq')
    
    xv = (-1+2*((0:N)+1/2)/(N+1));
    X = repmat(xv(:),1,N+1);
    Y = repmat(xv(:)',N+1,1);
    
    RXY = sqrt(X.^2+Y.^2);
    orig_img = 10*cos(3*RXY*pi/2);
    ind = find(RXY>1/2);
    orig_img(ind) = 10*cos(RXY(ind)*pi/2);
    ind = find(X>0 & X<.75 & Y>0 & Y<.75);
    orig_img(ind) = 10*sin(pi*X(ind)/2);
    imagesc(orig_img)
    
elseif strcmp(funct,'hill')
    
    xv = (-1+2*((0:N)+1/2)/(N+1));
    X = repmat(xv(:),1,N+1);
    Y = repmat(xv(:)',N+1,1);
    
    RXY = sqrt(X.^2+Y.^2);
    orig_img = 10*cos(3*RXY*pi/2);
    ind = find(RXY>1/2);
    orig_img(ind) = 10*cos(RXY(ind)*pi/2);
    
end

f = imresize(orig_img,[N,N]);

end

