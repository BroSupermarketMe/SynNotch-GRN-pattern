function X = padding(X,padsize)
    arguments
        X (1,:) struct
        padsize (1,1) double = 2
    end
    if mod(padsize,2)~=0
        error("padsize must be even")
    end

    % X is a 2D array: X.array0 and X.array1
    X.array0 = [zeros(size(X.array0,1),padsize),X.array0,zeros(size(X.array0,1),padsize)];
    X.array1 = [zeros(size(X.array1,1),padsize),X.array1,zeros(size(X.array1,1),padsize)];
    X.array0 = [zeros(padsize/2,size(X.array0,2));X.array0;zeros(padsize/2,size(X.array0,2))];
    X.array1 = [zeros(padsize/2,size(X.array1,2));X.array1;zeros(padsize/2,size(X.array1,2))];
    
end