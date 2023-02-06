function m = local_mass_1D(xcoord,nquad)
%
% LOCAL_MASS_1D - Function that returns local contribution of the
%                      1D advection matrix
%
% INPUT: 
%   xcoord - physical coordinates 
%   nquad  - number of quadrature points
%
% OUTPUT: 
%   m - local contribution of the mass matrix
%

%Get coordinates
x1 = xcoord(1);
x2 = xcoord(2);

%Size of the element
h = x2-x1;

%Quadrature points and weights
[z, w] = lobpts(nquad);

%Shape functions
N1 = @(x) (x2-x)/h;
N2 = @(x) (x-x1)/h;

%Jacobian
J = h/2;

%Initialize local matrix
m = zeros(2,2);

%Loop through quadrature points 
for i = 1:nquad
    
    %Get the evaluation point
    l = x1+(h/2)*(1+z(i));
    
    %Accumulate the local mass matrix   
    m = m+[N1(l)*N1(l) N1(l)*N2(l);...
         N2(l)*N1(l) N2(l)*N2(l)]*w(i)*abs(J);
end
end