% This script is to get Spherical Harmonics Descriptors
function run()

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
radii = output(:,1);
x1 = output(:,2);
y1 = output(:,3);
z1 = output(:,4);

%%
%GET SPHERICAL HARMONICS
% max radius is ceil(R*sqrt(2)) = 46; so we set m and r = 46
max_m = 46;
max_r = 46;
SH = zeros(max_r,max_m);

count = 1;
for iter = 1:max_r
    if(radii(count) < iter)
        count = count+1;
    end
end
count

end

%%
%PRE-PROCESSING
function [x_grid,y_grid,z_grid,final_dist] = pre_process(x1,y1,z1)
    %align coordinates
    x1 = x1-min(x1);
    y1 = y1-min(y1);
    z1 = z1-min(z1);

    %normalize to 1
    max_value = max([max(x1),max(y1),max(z1)]);
    x1 = x1./max_value;
    y1 = y1./max_value;
    z1 = z1./max_value;

    %rasterize to 2Rx2Rx2R voxel grid
    R = 32;
    x1 = round(x1*2*R);
    y1 = round(y1*2*R);
    z1 = round(z1*2*R);

    grid = zeros(2*R,2*R,2*R); 
    n_points = length(x1);
    for j = 1:n_points
        grid(y1(j)+1,x1(j)+1,z1(j)+1) = 1;
    end

    %get coordinates of grid voxel
    temp=reshape(grid,prod(size(grid)),1);
    [y_grid,x_grid,z_grid] = ind2sub(size(grid),find(temp==1));
    % x_grid = [];
    % y_grid = [];
    % z_grid = [];
    % for i = 1:2*R
    %     for j = 1:2*R
    %         for k = 1:2*R
    %             if(grid(i,j,k)==1)
    %                 y_grid = [y_grid; i];
    %                 x_grid = [x_grid; j];
    %                 z_grid = [z_grid; k];
    %             end
    %         end
    %     end
    % end

    %get CoM (center of mass)
    x_center = mean(x_grid);
    y_center = mean(y_grid);
    z_center = mean(z_grid);

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
    mean_dist = mean(dist);
    scale_ratio = (R/2)/mean_dist;
    x_grid = x_grid * scale_ratio;
    y_grid = y_grid * scale_ratio;
    z_grid = z_grid * scale_ratio;

    
    final_dist = sqrt((x_grid).^2 + (y_grid).^2 + (z_grid).^2);
    %verify
    % mean_dist = mean(dist)
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
