function [frames, pos, vel, velSummaries] = bird_flocks(numBirds, weights, params, viewSize, start, pred, iter, createMov)
%inputs: 
%         numBirds: number of birds to simulate
%          weights: rule weight vector; rules: separation, alignment, cohesion
%           params: parameter vector; params: r_s, v_m, r_p; r_p only needed if pred = 1
%        viewSize: size of viewing window
%           start: initialization setting
%                     1: start around left side, moving right, width will be twice viewSize
%                     2: start around center, moving upward
%            pred: whether to simulate a predator (1) or not (0)
%            iter: number of iterations to simulate
%       createMov: create a MATLAB movie
%
%outputs: 
%       frames: 3D matrix containing bird plot for each iteration, 1 indicates presence of a bird, 0 indicates empty space
%          pos: 2D matrix of positions; each row is one iteration, first two columns are x and y positions, in that order, of bird 1, second two are for bird 2 etc.
%          vel: 2D matrix of velocities, same structure as pos
% velSummaries: column vector showing average speed of birds at each iteration

r_s = params(1); % radius of separation
v_m = params(2); % maximum speed
r_p = params(3); % predator detection radius
m1 = weights(1);
m2 = weights(2);
m3 = weights(3);
if pred == 1
    m4 = weights(4);
end

frames = zeros(viewSize, viewSize, iter+1);
pos = zeros(iter+1, numBirds*2);
vel = zeros(iter+1, numBirds*2);
ruleCont = zeros(iter+1, 8);
velSummaries = zeros(iter+1, 1); 

predX = zeros(iter+1,2); % position of predator
predV = zeros(2,1); % velocity of predator


% initialize positions
if start == 1
    for i = 1:numBirds
        pos(1,2*i-1) = 0.2*viewSize*(rand() + 1);  %x positions
        vel(1,2*i-1) = 0.04*viewSize;
        pos(1,2*i) = floor(viewSize/2) + 0.1*viewSize*rand(); %y positions
        vel(1,2*i) = 0.02*viewSize*(2*rand()-1);
        if pred == 1
            predX(1,:) = [viewSize/2, 0];
            predV = [0.1*viewSize, 0.1*viewSize];
        end
    end
elseif start == 2
    for i = 1:numBirds
        pos(1,2*i-1) = viewSize/2 + 0.1*viewSize*(2*rand() - 1);  %x positions
        vel(1,2*i-1) = 0;
        pos(1,2*i) = viewSize/2 + 0.1*viewSize*(2*rand() - 1); %y positions
        vel(1,2*i) = 0.05*viewSize*(2*rand()-1);
        if pred == 1
            predX(1,:) = [0,0];
            predV = [0.075*viewSize, 0.075*viewSize];
        end
    end
end

% simulate
for i = 2:iter+1
    for j = 1:numBirds
        if pred == 1
            if rule4(j,i-1) > 0
                delta_v = m4.*rule4(j,i-1);
            else
                delta_v = m1.*rule1(j,i-1) + m2.*rule2(j,i-1) + m3.*rule3(j,i-1);
            end
        else 
            delta_v = m1.*rule1(j,i-1) + m2.*rule2(j,i-1) + m3.*rule3(j,i-1);
        end
        if j == 10
            ruleCont(i, 1:2) = rule1(j,i-1);
            ruleCont(i, 3:4) = rule2(j,i-1);
            ruleCont(i, 5:6) = rule3(j,i-1);
            ruleCont(i, 7:8) = rule4(j,i-1);
        end
        v = capVecMag([vel(i-1,2*j-1) + delta_v(1), vel(i-1,2*j) + delta_v(2)], v_m);
        vel(i,2*j-1) = v(1);
        vel(i,2*j) = v(2);
        pos(i,2*j-1) = pos(i-1,2*j-1) + vel(i,2*j-1); 
        pos(i,2*j) = pos(i-1,2*j) + vel(i,2*j); 
    end
    predX(i,:) = predX(i-1,:) + predV;

    k = arrayfun(@(x,y) sqrt(x^2 + y^2), vel(i,1:2:end), vel(i,2:2:end));
    velSummaries(i,1) = mean(k);
end

% generate frames
if createMov == 1
    for i = 1:iter+1
        for j = 1:numBirds
            x = pos(i,2*j-1);
            y = pos(i,2*j);
            prX = predX(i,1);
            prY = predX(i,2);
            [r, c] = coordsToIndices(prX, prY);
            frames(r,c,i) = 1;
            [r,c] = coordsToIndices(x,y);
            frames(r,c,i) = 1;
        end
    end

    M(iter+1) = struct('cdata',[],'colormap',[]);
    for i = 1:iter+1
        [r,c] = find(frames(:,:,i));
        figure
        scatter(c, r, 75, 'sw', 'filled')
        set(gca, 'Color','k', 'YDir','reverse')
        if start == 1
            k = 2;
        elseif start == 2
            k = 1;
        end
        axis([0 k*viewSize  0 viewSize])
        drawnow
        M(i) = getframe;
    end
    movie(M);
    r = int2str(randi([1,1000], 1,1));
    name = strcat('sim_',r,'.avi');
    v = VideoWriter(name, 'Motion JPEG AVI');
    open(v);
    writeVideo(v, M);
    close(v);
    close all
    disp(strcat('movie is',{' '}, name))
end

disp(ruleCont);



% separation
function k1  = rule1(i, k) % i is bird index, k is iteration to evaluate
    k1 = zeros(1,2);
    x = pos(k, 2*i-1); % position of bird i
    y = pos(k, 2*i);
    for a = 1:numBirds
        if i == a
            continue
        end
        thisX = pos(k,2*a-1); % position of bird a
        thisY = pos(k,2*a);
        dist = sqrt((x-thisX)^2 + (y-thisY)^2);
        if dist < r_s
            vec = [x,y] - [thisX,thisY];
            dir = vec/norm(vec);
            k1 = k1 + (1 - dist/r_s).*dir;
        end
    end
end

% alignment
function k2 = rule2(i, k) % i is bird index, k is iteration to evaluate
    sum = zeros(1,2);
    for a = 1:numBirds
        if i == a
            continue
        end
        sum = sum + [vel(k, 2*a-1) vel(k, 2*a)];
    end
    k2 = sum / (numBirds - 1);
end

% cohesion
function k3 = rule3(i, k) % i is bird index, k is iteration to evaluate
    sum = zeros(1,2);
    for a = 1:numBirds
        if i == a
            continue
        end
        sum = sum + [pos(k, 2*a-1) pos(k, 2*a)];
    end
    com = sum / (numBirds - 1);
    k3 = com - [pos(k,2*a-1) pos(k,2*a)];
end

% predator avoidance
function k4 = rule4(i, k) % i is bird index, k is iteration to evaluate
    k4 = zeros(1,2);
    x = pos(k, 2*i-1); % position of bird i
    y = pos(k, 2*i);
    prX = predX(k,1); % position of predator
    prY = predX(k,2);
    dist = sqrt((x-prX)^2 + (y-prY)^2);
    if dist < r_p
        vec = [x,y] - [prX, prY];
        dir = vec/norm(vec);
        k4 = (1 - dist/r_p).*dir;
    end
end


function [r,c] = coordsToIndices(x,y)
    if y < 0
        y = 0;
    elseif y > viewSize-1
        y = viewSize - 1;
    end
    r = round(viewSize - y);
    if x < 1
        x = 1;
    elseif x > 2*viewSize
        x = 2*viewSize;
    end
    c = round(x);
end

function [vsc] = capVecMag(v, lim) %  if the magnitude of v is greater than lim, rescale v so it has magnitude lim
    if norm(v) > lim
        vsc = lim.*(v / norm(v));
    else 
        vsc = v;
    end
end

end