function g=local_advection_1D(xcoord,beta,nquad)
%
% LOCAL_ADVECTION_1D - Function that returns local contribution of the
%                      1D advection matrix
%
% INPUT: 
%   xcoord - physical coordinates 
%   beta   - 1D component of the advection vector
%   nquad  - number of quadrature points
%
% OUTPUT: 
%   g - local contribution of the advection matrix
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

%Derivatives of shape functions
dN1 = @(x) -1/h;
dN2 = @(x) 1/h;

%Jacobian
J = h/2;

%Initialize local matrix
g = zeros(2,2);

%Loop through quadrature points 
for i = 1:nquad
    
    %Get the evaluation point
    l = x1+(h/2)*(1+z(i));
    
    %Accumulate the local advection matrix
    g = g+[dN1(l)*N1(l) dN1(l)*N2(l);...
         dN2(l)*N1(l) dN2(l)*N2(l)]*beta*w(i)*abs(J);
end

end


