function [a,r,c]=ASAtransform(p0,q0,coord_type)
    arguments
        p0 (1,1) double
        q0 (1,1) double
        coord_type string = "pq"
    end
    if coord_type == "pq"
        c=p0-1;
        a=mod(q0-1,2);
        r=floor((q0-1)/2);
    elseif coord_type == "xy"
        c=q0;
        a=mod(p0,2);
        r=floor(p0/2);
    end
end