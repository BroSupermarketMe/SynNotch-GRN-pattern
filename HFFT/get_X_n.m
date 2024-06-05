function [X_n] = get_X_n(X)
    r0=size(X.array0,1);
    r1=size(X.array1,1);
    if not ((r0==r1) & (mod(r0,2)==0))
        error('X.array0 and X.array1 must have the same even number of rows');
    end

    % 获得F的每个元素的下标
    [m,n] = size(X.array0);
    X_n = zeros(2*m,n);
    for i = 1:m
        for j = 1:n
            % 以复数的形式存储
            [x,y] = get_x_y(0,i,j);
            X_n(2*i-1,j) = complex(x,y);
            [x,y] = get_x_y(1,i,j);
            X_n(2*i,j) = complex(x,y);
        end
    end
    [x_center,y_center] = get_x_y(0,m/2+1,n/2+1);
    % 将位置矩阵居中
    X_n = X_n - complex(x_center,0);
    X_n = X_n - complex(0,y_center);
end

function [x,y]=get_x_y(a,r,c)
    T_matrix=[1/2 0 1; sqrt(3)/2 sqrt(3) 0];
    xy=T_matrix*[a;r;c];
    x=xy(1);
    y=xy(2);
end