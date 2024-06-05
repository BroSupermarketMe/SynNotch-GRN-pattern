function X=hfftshift(X)
    % HFFTSHIFT Shifts the zero frequency component to the center of the spectrum
    % X: structure with fields: array0, array1
    r0=size(X.array0,1);
    r1=size(X.array1,1);
    if r0~=r1
        error('The two arrays must have the same size to be shifted together.');
    end
    if mod(r0,2)~=0
        error('The size of the arrays must be even.');
    end
    X.array0=fftshift(X.array0);
    X.array1=fftshift(X.array1);
end