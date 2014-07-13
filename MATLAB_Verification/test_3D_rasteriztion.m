%rasterization between two 3D points
function run()
    clc
    clear
    
    grid = zeros(10,10,10);
    p1 = [1,1,1];
    p2 = [10,10,10];
    p3 = [1,10,10];
    [len12,coord12,grid] = fillGridLine(p1,p2,grid);
    [len13,coord13,grid] = fillGridLine(p1,p3,grid);
    [len23,coord23,grid] = fillGridLine(p2,p3,grid);
    
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
    x_grid
    y_grid
    z_grid
    
    %draw
   scatter3(x_grid,y_grid,z_grid,50,[0 0 1],'*'); view([60,-60,60]);
    
    
end

function [len,coord,grid] = fillGridLine(p1,p2,grid)

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