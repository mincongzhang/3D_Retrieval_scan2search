%rasterization between two 3D points
function run()
    clc
    clear
    close all
    grid = zeros(10,10,10);
    p1 = [1,1,1];
    p2 = [10,10,10];
    p3 = [1,10,10];

    [grid,coord12,len12] = fillGridLine(p1,p2,grid);
    [grid,coord23,len23] = fillGridLine(p2,p3,grid);
    [grid,coord13,len13] = fillGridLine(p1,p3,grid);

    len12
    len13
    diff = len12-len13;
    if(diff>=0)
        for i = 1:len13
            i
            [grid,~,~] = fillGridLine(coord12(i,:),coord13(i,:),grid);
        end
        
        for i = len13+1:len12
            i
            [grid,~,~] = fillGridLine(coord12(i,:),coord13(len13,:),grid);
        end
        
    end
    
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

function [grid,coord,len] = fillGridLine(p1,p2,grid)

    v12 = p2-p1;
    len_p12 = sqrt(v12(1)^2+v12(2)^2+v12(3)^2);
    v12_norm = v12./(len_p12);
    
    coord = [];
    len = 0;
    
    for i = 0:floor(len_p12)
        tmp_p = round(p1+v12_norm*i);
        if(grid(tmp_p(1),tmp_p(2),tmp_p(3))~=1)
            grid(tmp_p(1),tmp_p(2),tmp_p(3))=1;
            coord = [coord;tmp_p(1),tmp_p(2),tmp_p(3)];
            len = len+1;
        end
    end
    grid(p2(1),p2(2),p2(3))=1;
end