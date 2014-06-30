% This script is to get Spherical Harmonics Descriptors in a regular way
function F_r = SH()

close all
fv = stlread('femur.stl');
vertices = fv.vertices;

x1 = vertices(:,1);
y1 = vertices(:,2);
z1 = vertices(:,3);

%%
%PRE-PROCESSING
[x_grid,y_grid,z_grid,dist] = pre_process(x1,y1,z1);

%%
%SORT ACCORDING TO DISTANCE/RADIUS (MIN to MAX)
% size(x_grid) = (9690,1)
output = qsort_radius(dist',x_grid',y_grid',z_grid');
output = output';

% size(output) = (9690,4)
dist_vector = output(:,1);
x     = output(:,2);
y     = output(:,3);
z     = output(:,4);
phi   = atan(y./x);
theta = acos(z./dist_vector);
phi(1)
theta(1)

%test SH for same theta and phi, with different radius
% radii = [1:32]';
% phi = [ones(16,1)*pi/2;ones(16,1)*pi];
% theta = [ones(16,1)*pi/2;ones(16,1)*pi];

%%
%GET SPHERICAL HARMONICS
max_l = 32;
max_r = 32;
F_r = zeros(max_r,1);
SH = zeros(max_r,max_l);

idx_n = 1;
for idx_r = 1:max_r
    while(dist_vector(idx_n) < idx_r) %BUGGGGGGGGGGG HERE???
        % calculate SH
        F_lr = zeros(max_l,1);
        
        %for each F_r
        for idx_l = 1:max_l
            Y_ml = zeros(2*idx_l+1,1);
            ml_count = 1;
            %for each F_lr
            for idx_m = -idx_l:idx_l
                if(idx_m>=0)
                    Y_ml(ml_count,1) = spharm(idx_l,idx_m,theta(idx_n),phi(idx_n));
                else
                    Y_temp = spharm(idx_l,-idx_m,theta(idx_n),phi(idx_n));
                    Y_ml(ml_count,1) = (-1)^(-idx_m) * conj(Y_temp);
                end
                ml_count = ml_count+1;
            end
            F_lr(idx_l,1) = sum(Y_ml);
        end
        F_r(idx_r,1) = sum(F_lr);
        sum(F_lr)
        
        if(idx_n==length(dist_vector))
            break;
        else
            idx_n = idx_n+1
        end
    end
end

end

%%
%PRE-PROCESSING
function [x_grid,y_grid,z_grid,final_dist] = pre_process(x1,y1,z1)
    %align coordinates
    x1 = x1-min(x1);
    y1 = y1-min(y1);
    z1 = z1-min(z1);

    %normalize to 1 and rasterize to 2Rx2Rx2R voxel grid
    max_value = max([max(x1),max(y1),max(z1)]);
    R = 32;
    x1 = round(x1./max_value*(2*R-1));
    y1 = round(y1./max_value*(2*R-1));
    z1 = round(z1./max_value*(2*R-1));
    
    x_grid_grid = zeros(2*R);
    y_grid_grid = zeros(2*R);
    z_grid_grid = zeros(2*R);
    x_grid = [];
    y_grid = [];
    z_grid = [];
    grid = zeros(2*R,2*R,2*R);
    n_points = length(x1)
    for j = 1:n_points
        if(grid(x1(j)+1,y1(j)+1,z1(j)+1)==0)
            %register
            grid(x1(j)+1, y1(j)+1, z1(j)+1)=1;
            x_grid = [x_grid x1(j)];
            y_grid = [y_grid y1(j)];
            z_grid = [z_grid z1(j)];
        end
    end   
    x_grid = x_grid';
    y_grid = y_grid';
    z_grid = z_grid';
    x_grid_size = size(x_grid)

    %get CoM (center of mass)
    x_center = mean(x_grid)
    y_center = mean(y_grid)
    z_center = mean(z_grid)

    %move CoM to (0,0,0)
    x_grid = x_grid - x_center;
    y_grid = y_grid - y_center;
    z_grid = z_grid - z_center;

    %verify
    % mean(x_grid)
    % mean(y_grid)
    % mean(z_grid)

    %scale and make the average distance to CoM is R/2;
    dist = sqrt((x_grid).^2 + (y_grid).^2 + (z_grid).^2);
    mean_dist = mean(dist)
    scale_ratio = (R/2)/mean_dist
    x_grid = x_grid * scale_ratio;
    y_grid = y_grid * scale_ratio;
    z_grid = z_grid * scale_ratio;

    
    final_dist = sqrt((x_grid).^2 + (y_grid).^2 + (z_grid).^2);
    %verify
    final_mean_dist = mean(final_dist)
    % 
    % scatter3(x1,y1,z1,5,[0 0 1],'.'); view([60,-60,60]);
    % set(gca, 'XLim', [-100 100]);
    % set(gca, 'YLim', [-100 100]);
    % set(gca, 'ZLim', [-100 100]);
    % figure,
    % scatter3(x_grid,y_grid,z_grid,5,[0 0 1],'.'); view([60,-60,60]);
    % set(gca, 'XLim', [-100 100]);
    % set(gca, 'YLim', [-100 100]);
    % set(gca, 'ZLim', [-100 100]);
end

%%
%SORT ACCORDING TO DISTANCE/RADIUS (MIN to MAX)
function output = qsort_radius(list,x,y,z)
    startpoint = 1;
    endpoint = length(list);
    output = [list;x;y;z];
    if(startpoint < endpoint)

        flag = startpoint;
        for j = (startpoint+1):endpoint
            if(output(1,startpoint)>output(1,j))
                flag = flag+1;
                %swap
                temp = output(:,flag);
                output(:,flag) = output(:,j);
                output(:,j) = temp;
            end
        end
        %swap
        temp = output(:,startpoint);
        output(:,startpoint) = output(:,flag);
        output(:,flag) = temp;
        
        %recursive
        output(:,startpoint:flag) = qsort_radius(output(1,startpoint:flag),output(2,startpoint:flag),output(3,startpoint:flag),output(4,startpoint:flag));
        output(:,flag+1:endpoint) = qsort_radius(output(1,flag+1:endpoint),output(2,flag+1:endpoint),output(3,flag+1:endpoint),output(4,flag+1:endpoint));
    end
end

%%
%GET SPHERICAL HARMONICS
