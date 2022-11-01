function [bn,berr,bl,bcost] = setup_quadrature(normA, tol, maxp, method)
%SETUP_QUADRATURE    Estimates the optimal scaling and number of quadrature nodes for phiquadmv
%  [bn,berr,bl,bcost] = setup_quadrature(normA)
%  computes the optimal scaling factor bl, number of quadrature nodes bn, and related a priori
%  error bound berr and total cost bcost that minimizes the total cost while respecting
%  a default quadrature error tolerance of 10^{-16} for the Gauss-Legendre version of phiquadmv
%  with q <= 20. The input normA must be an estimate of the norm of A.
%
%  [bn,berr,bl,bcost] = setup_quadrature(normA, tol)
%  is the same, but with a different tolerance tol.
%  
%  [bn,berr,bl,bcost] = setup_quadrature(normA, tol, maxp)
%  is the same, but for all q <= maxp
%
%  [bn,berr,bl,bcost] = setup_quadrature(normA, tol, maxp, 'gauss')
%  is the same as the above.
%
%  [bn,berr,bl,bcost] = setup_quadrature(normA, tol, maxp, method)
%  for any string different than 'gauss' is the same, but for Clenshaw-Curtis quadrature.

    if nargin < 4
        method = 'gauss';
    end
    if nargin < 3
        maxp = 20;
    end
    if nargin < 2
        tol = 1e-16;
    end

    maxp = maxp + 1;

    minl = 0;
    maxl = ceil(log2(normA));
    get_cost = @(n,l) n + l*maxp;

    bn = inf; berr = inf; bl = inf; bcost = inf;
    for l = maxl:-1:minl
        [n,err] = estimate_n(2^(-l)*normA, tol, maxp, method);
        cost = get_cost(n,l);
        if cost < bcost
            bn = n+1; berr = err; bl = l; bcost = cost;
        else
            break;
        end
    end
end

function [n,err] = estimate_n(normA, tol, maxp, method)
    if nargin < 4
        method = 'gauss';
    end
    if nargin < 3
        maxp = 20;
    end
    if nargin < 2
        tol = 1e-16;
    end

    logtol = log(tol);
    n = 1;
    err = get_errlog(n, normA, maxp, method);

    while err > logtol
        n = 2*n;
        err = get_errlog(n, normA, maxp, method);
    end

    if n > 1
        [res,fval,exitflag,output] = fzero(@(n) get_errlog(n, normA, maxp, method)-logtol, [n/2, n]);
        if exitflag > 0
            n = ceil(res);
            err = get_errlog(n, normA, maxp, method);
        end
    end
end

function err = get_errlog(n, normA, maxp, method)
    if nargin < 3
        maxp = 20;
    end
    if nargin < 4
        method = 'gauss';
    end

    err = -inf;
    for p = 1:maxp
        if strcmp(method, 'gauss')
            rs = roots([normA/4, p-2-2*n, -normA/2-2*p, 2*n + p, normA/4]);
            rs = real(rs(isreal(rs) & real(rs) > 1));
            temp = min(errlog(rs, n, p, normA));
        else
            rs = roots([normA/4, p-1-n, -normA/2-2*p, n-1 + p, normA/4]);
            rs = real(rs(isreal(rs) & real(rs) > 1));
            temp = min(errlog_cc(rs, n, p, normA));
        end
        err = max(err, temp);
    end
end

function err=errlog(r, n, p, normA)
    g = @(r) (r+1).^2./(2*r);
    err = log(144/35) - log(r.^2-1) - 2*n*log(r) + p*log(g(r))-(p+1)*log(2)-gammaln(p+1)+0.5*g(r)*normA;
end

function err=errlog_cc(r, n, p, normA)
    g = @(r) (r+1).^2./(2*r);
    err = log(144/35) - log(r.^2-1) + (1-n)*log(r) + p*log(g(r))-(p+1)*log(2)-gammaln(p+1)+0.5*g(r)*normA;
end
