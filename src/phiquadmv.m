function [y,varargout] = phiquadmv(q,As,b,normA,scaling,varargin)
%PHIQUADMV   Computes the phi matrix functions appearing in expint methods.
%   y = phiquadmv(q, As, b, normA, scaling)
%   or y = phiquadmv(q, As, b, normA, scaling, 'chebadaptive')
%   returns phi_q(A)*b for a vector b, where A
%   has Kronecker sum structure and is expressed in terms of 1D matrices
%   contained in a length 2 (in 2D) or length 3 (in 3D) cell array As, i.e.
%   A = As{1} \oplus As{2}, (in 2D), or A = As{1} \oplus As{2} \oplus As{3} (in 3D).
%   If As is prescribed as a single 1D matrix then it is assumed that
%   A = \bigoplus_{i=1}^d As, where d is the dimension which is computed by checking
%   the size of b. If q > 1, then the returned y is of size lenght(b)-by-q containing
%   phi_i for i=1,...,q. normA must be a scalar being an estimate of norm(A,1). 
%   scaling must be an integer so that the matrix will be scaled by 2^{scaling}.
%   If the scaling prescribed is negative, then the optimal scaling factor will be computed
%   given a default error tolerance of 10^{-14}. y is computed by using adaptive
%   Clenshaw-Curtis quadrature by default.
%
%   y = phiquadmv(q, As, b, normA, scaling, tol) for scaling < 0 uses Gauss-Legendre quadrature
%   instead with the scaling and number of points chosen optimally to satisfy an error tolerance of
%   tol and minimize the total cost. Gauss-Legendre quadrature is slightly more efficient so this option
%   is to be preferred.
%
%   y = phiquadmv(q, As, b, normA, scaling, npts) for scaling >= 0 and npts integer uses Gauss-Legendre quadrature
%   instead with the prescribed scaling and number of quadrature points.
%
%   y = phiquadmv(q, As, b, normA, scaling, 'gauss', tol) and y = phiquad(q, As, b, normA, scaling, 'gauss', npts)
%   are the same as the two above depending on the value of scaling.
%
%   y = phiquadmv(q, As, b, normA, scaling, 'chebadaptive', tol) uses adaptive Clenshaw-Curtis with relative error
%   tolerance tol and scaling factor given by scaling if this is non-negative. If scaling < 0, then the optimal
%   scaling factor will be estimated.
%
%   Additional outputs are:
%   [y,npts,err,nkron,l] = phiquadmv(q, As, b, normA, scaling),
%   where npts is the final number of quadrature nodes used, 
%   nkron is the number of actions of exp(A) computed, l is the final scaling used,
%   and err is the final estimated relative error
%
%   ||y_{2*npts+1}-y_{npts+1}||/||y_{npts+1}||
%
%   where y_{m} is the result computed using m quadrature nodes (Clenshaw-Curtis only, otherwise nan is returned).
%   Any number of outputs between 0 and 5 can be returned.
%   
%   NOTE: requires chebfun to be installed and in the path.
%         See https://www.chebfun.org/
%
%   Copyright @ 2022 by Matteo Croci and Judit Muñoz-Matute.

    norb = norm(b,inf);
    b = b/norb;

    tol = 1.0e-14;
    nout = max(nargout,1) - 1;
    err = nan;

    if iscell(As)
        dim = length(As);
    else
        n = size(As,1);
        nb = length(b);
        if n^2 == nb
            dim = 2;
        elseif n^3 == nb
            dim = 3;
        else
            error('Only implemented for dimensions 2 and 3.');
        end
    end
    
    if dim == 2
        expmatvec_fun = @(s,b) expmatvec2D(s,b,As);
    elseif dim == 3
        expmatvec_fun = @(s,b) expmatvec3D(s,b,As);
    else
        error('Only implemented for dimensions 2 and 3.');
    end

    if q == 0
        y = norb*expmatvec_fun(1,b);
        return
    else
        q = 1:max(q);
    end

    lq = length(q);

    if nout > 4
        error('Too many outputs requested!')
    end

    if isempty(varargin) || (length(varargin) == 1 && isstring(varargin{1}))
        method = 'chebadaptive';
    elseif length(varargin) == 1 && ~isstring(varargin{1})
        method = 'gauss';
        if scaling >= 0
            npts = varargin{1};
        else
            tol = varargin{1};
        end
    elseif length(varargin) == 2
        method = varargin{1};
        if strcmp(method, 'gauss') 
            if scaling >= 0
                npts = varargin{2};
            else
                tol = varargin{2};
            end
        else
            tol = varargin{2};
        end
    else
        error('Invalid arguments!')
    end

    if scaling >= 0
        l = round(scaling);
    elseif strcmp(method, 'gauss')
        [npts,~,l,~] = setup_quadrature(normA, tol/norb, max(q), method);
    else 
        [~,~,l,~] = setup_quadrature(normA, tol/norb, max(q), method);
    end

    scaled_fun = @(s) expmatvec_fun(s*2^(-l),b);

    if ~strcmp(method, 'chebadaptive')
        [x,w] = legpts(npts,[0,1]);

        V = zeros(length(b), npts);
        for i=1:npts-1
            V(:,i) = scaled_fun(1-x(i));
        end

        V(:,end) = scaled_fun(1-x(end));

        y = zeros(length(b),lq);
        for j=1:lq
            y(:,j) = V*(w'.*x.^(q(j)-1)/factorial(q(j)-1));
        end

    else
        n = 4;
        [x,w] = chebpts(2*n+1,[0,1]);
        
        V = zeros(length(b),2*n+1);
        for i=1:2*n
            V(:,i) = scaled_fun(1-x(i));
        end
        V(:,end) = b;
        
        err = inf;
        while err > tol && n < 1000

            [x0,w0] = chebpts(n+1,[0,1]);

            %if ~isnested(x,x0)
            %    error('not nested')
            %end

            y0 = zeros(length(b),lq);
            y = zeros(length(b),lq);

            for j=1:lq
                y(:,j) = V*(w'.*x.^(q(j)-1)/factorial(q(j)-1));
                y0(:,j) = V(:,1:2:2*n+1)*(w0'.*x0.^(q(j)-1)/factorial(q(j)-1));
            end

            err = max(max(abs(y-y0)))/max(max(abs(y)));

            n = 2*n;
            
            [x,w] = chebpts(2*n+1,[0,1]);

            VV = zeros(length(b),2*n+1);
            VV(:,1:2:2*n+1) = V;
            for i=2:2:2*n
                VV(:,i) = scaled_fun(1-x(i));
            end
            V = VV;

        end

        npts = 2*n+1;
    end
    
    if l > 0
        y = norb*square_back(y,l, As);
    end

    if nout == 1
        varargout{1} = npts;
    elseif nout == 2
        varargout{1} = npts;
        varargout{2} = err;
    elseif nout == 3
        varargout{1} = npts;
        varargout{2} = err;
        varargout{3} = npts + l*max(q);
    elseif nout == 4
        varargout{1} = npts;
        varargout{2} = err;
        varargout{3} = npts + l*max(q);
        varargout{4} = l;
    end

end

function ys = square_back(y, l, As)
    % modified squaring algorithm
    [nb,lq] = size(y);
    isc = iscell(As);

    if isc
        dim = length(As);
    else
        n = size(As,1);
        if n^2 == nb
            dim = 2;
        elseif n^3 == nb
            dim = 3;
        else
            error('Only implemented for dimensions 2 and 3.');
        end
    end

    if isc
        if dim == 2
            eAs = {expm(2^(-l)*As{1}), expm(2^(-l)*As{2})};
            readyexpmatvec = @(eAs,b) readyexpmatvec2D(eAs, b);
        else
            eAs = {expm(2^(-l)*As{1}), expm(2^(-l)*As{2}), expm(2^(-l)*As{3})};
            readyexpmatvec = @(eAs,b) readyexpmatvec3D(eAs, b);
        end
    else
        eAs = expm(2^(-l)*As);
        if dim == 2
            readyexpmatvec = @(eAs,b) readyexpmatvec2D(eAs, b);
        else
            readyexpmatvec = @(eAs,b) readyexpmatvec3D(eAs, b);
        end
    end

    ys = y;
    for it=1:l
        for i=1:lq
            ys(:,i) = ( readyexpmatvec(eAs,y(:,i)) + y(:,1:i)*(1./factorial(i-(1:i)')) )/2^i;
        end
        if isc
            for i=1:dim
                eAs{i} = eAs{i}*eAs{i};
            end
        else
            eAs = eAs*eAs;
        end
        y = ys;
    end
end

function y = readyexpmatvec2D(eAs, b)
    if iscell(eAs)
        m1 = size(eAs{1},1);
        m2 = size(eAs{2},1);
        y = reshape(eAs{1}*reshape(b,m1,m2)*eAs{2}',m1*m2,1);
    else
        m = size(eAs,1);
        y = reshape(eAs*reshape(b,m,m)*eAs',m*m,1);
    end
end

function y = expmatvec2D(s,b,As)
    if iscell(As)
        eAs = {expm(s*As{1}), expm(s*As{2})};
        y = readyexpmatvec2D(eAs,b);
    else
        eA = expm(s*As);
        y = readyexpmatvec2D(eA,b);
    end
end

function y = readyexpmatvec3D(eAs, b)
    if iscell(eAs)
        m1 = size(eAs{1},1);
        m2 = size(eAs{2},1);
        m3 = size(eAs{3},1);
        y = reshape(tucker(reshape(b,m1,m2,m3), eAs{1}, eAs{2}, eAs{3}), m1*m2*m3, 1);
    else
        m = size(eAs,1);
        y = reshape(tucker(reshape(b,m,m,m), eAs, eAs, eAs), m*m*m, 1);
    end
end

function y = expmatvec3D(s,b,As)
    if iscell(As)
        eAs = {expm(s*As{1}), expm(s*As{2}), expm(s*As{3})};
        y = readyexpmatvec3D(eAs,b);
    else
        eA = expm(s*As);
        y = readyexpmatvec3D(eA,b);
    end
end
