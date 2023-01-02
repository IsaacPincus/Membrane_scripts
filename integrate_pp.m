function int = integrate_pp(pp,breaks)
    % integrates a pp form cubic spline, which is the output from csape()
    
    if nargin == 1
        breaks = pp.breaks;
    end

    n = pp.pieces;
    int = zeros(1,n+1);

    for ii=1:n
        min_x = breaks(ii);
        max_x = breaks(ii+1);
        p = pp.coefs(ii,:);
        q = polyint(p);
        int(ii+1) = int(ii) + diff(polyval(q,[0, max_x-min_x]));
    end

end