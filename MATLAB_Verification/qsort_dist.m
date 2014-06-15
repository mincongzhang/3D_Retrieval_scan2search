function run()
clc
list = [4 7 5 3 9 1 4 6 3 5];
x = [1 2 3 4 5 6 7 8 9 10];
y = [1 2 3 4 5 6 7 8 9 10];
z = [1 2 3 4 5 6 7 8 9 10];
startpoint = 1;
endpoint   = length(list);
output = qsorting(list,x,y,z)

end

function output = qsorting(list,x,y,z)
    startpoint = 1;
    endpoint = length(list);
    output = [list;x;y;z];
    if(startpoint < endpoint)

        flag = startpoint;
        for i = (startpoint+1):endpoint
            if(output(1,startpoint)>output(1,i))
                flag = flag+1;
                %swap
                temp = output(:,flag);
                output(:,flag) = output(:,i);
                output(:,i) = temp;
            end
        end
        %swap
        temp = output(:,startpoint);
        output(:,startpoint) = output(:,flag);
        output(:,flag) = temp;
        
        %recursive
        output(:,startpoint:flag) = qsorting(output(1,startpoint:flag),output(2,startpoint:flag),output(3,startpoint:flag),output(4,startpoint:flag));
        output(:,flag+1:endpoint) = qsorting(output(1,flag+1:endpoint),output(2,flag+1:endpoint),output(3,flag+1:endpoint),output(4,flag+1:endpoint));
    end
end