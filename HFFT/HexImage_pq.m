function X = HexImage_pq(yout,P,Q)
    arguments
        yout (1,:) double
        P (1,1) double = 10
        Q (1,1) double = 10
    end
    [a,r_max,c_max] = ASAtransform(P,Q);
    if a == 0
        X.array0=zeros(r_max+1,c_max+1);
        X.array1=zeros(r_max,c_max+1);
    elseif a == 1
        X.array0=zeros(r_max+1,c_max+1);
        X.array1=zeros(r_max+1,c_max+1);
    else
        error("Invalid ASAtransform")
    end
    for i=1:length(yout)
        [p,q]=ind2pq(i,P);
        [a,r,c]=ASAtransform(p,q);
        if a == 0
            X.array0(r+1,c+1)=yout(i);
        elseif a == 1
            X.array1(r+1,c+1)=yout(i);
        else
            error("Invalid ASAtransform")
        end
    end
end

function [p,q]=ind2pq(ind, P)
    q = 1+floor((ind-1)/P);
    p = ind - (q-1)*P;
end