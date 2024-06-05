function yout = PQImage_hex(X)
    [r0_max,c0_max]=size(X.array0);
    [r1_max,c1_max]=size(X.array1);
    P=c0_max;
    if r0_max < r1_max
        Q=r0_max*2-1;
    else
        Q=r0_max*2;
    end
    yout = zeros(1,P*Q);
    for i=1:r0_max
        for j=1:c0_max
            [p,q]=PQtransform(0,i-1,j-1);
            ind=pq2ind(p,q,P);
            yout(ind)=X.array0(i,j);
        end
    end
    for i=1:r1_max
        for j=1:c1_max
            [p,q]=PQtransform(1,i-1,j-1);
            ind=pq2ind(p,q,P);
            yout(ind)=X.array1(i,j);
        end
    end
end

function ind=pq2ind(p,q, P)
    ind = p + (q-1)*P;
end