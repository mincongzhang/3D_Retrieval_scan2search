%rasterization between two 3D points
function run()
    clc
    clear
    
    grid = zeros(10,10,10);
    p1 = [1,1,1];
    p2 = [9,9,9];
    p3 = [1,9,10];
    max_x = max([p1(1),p2(1),p3(1)]);
    max_y = max([p1(2),p2(2),p3(2)]);
    max_z = max([p1(3),p2(3),p3(3)]);
    min_x = min([p1(1),p2(1),p3(1)]);
    min_y = min([p1(2),p2(2),p3(2)]);
    min_z = min([p1(3),p2(3),p3(3)]);
    x_len = max_x - min_x
    y_len = max_y - min_y
    z_len = max_z - min_z
    
    %返回最大len所在轴
    max_on_axis = 1;
    max_len = x_len;
    if(y_len>max_len)
        max_on_axis = 2;
        max_len = y_len; 
    end
    if(z_len>max_len)
        max_on_axis = 3;
        max_len = z_len; 
    end    
    max_len
    max_on_axis
    
    %排序，返回按其轴从大到小排列的3个点(C++里用stl的sort)
    tmp = [p1(max_on_axis);p2(max_on_axis);p3(max_on_axis)];
    pt_tmp = [p1;p2;p3]
    pt = [];
    for i = 1:3
       idx = find(tmp == max(tmp))
       tmp(idx)=[];
       pt = [pt;pt_tmp(idx,:)]; 
    end
    pt
    
    
    [coord12,grid] = fillGridLine(pt(1,:),pt(2,:),grid);
    [coord13,grid] = fillGridLine(pt(2,:),pt(3,:),grid);
    [coord23,grid] = fillGridLine(pt(1,:),pt(3,:),grid);
    coord12
    coord13
    coord23
    size(coord12)
    size(coord13)
    size(coord23)
    
    %get coordinates
    x_grid = [];
    y_grid = [];
    z_grid = [];
    for x = 1:10
        for y = 1:10
            for z = 1:10
                if(grid(x,y,z)==1)
                    x_grid = [x_grid x];
                    y_grid = [y_grid y];
                    z_grid = [z_grid z];
                end
            end
        end
    end
%     x_grid
%     y_grid
%     z_grid
    
    %draw
   scatter3(x_grid,y_grid,z_grid,50,[0 0 1],'*'); view([60,-60,60]);
    
    
end

function [coord,grid] = fillGridLine(p1,p2,grid)

    v12 = p2-p1;
    len_p12 = sqrt(v12(1)^2+v12(2)^2+v12(3)^2);
    v12_norm = v12./(len_p12);
    
    coord = [];
    for i = 0:floor(len_p12)
        tmp_p = round(p1+v12_norm*i);
        if(grid(tmp_p(1),tmp_p(2),tmp_p(3))~=1)
            grid(tmp_p(1),tmp_p(2),tmp_p(3))=1;
            coord = [coord;tmp_p(1),tmp_p(2),tmp_p(3)];
        end
    end
    grid(p2(1),p2(2),p2(3))=1;
end