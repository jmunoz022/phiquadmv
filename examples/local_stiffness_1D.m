function k = local_stiffness_1D(xcoord,diff,nquad)
%
% LOCAL_STIFFNESS_1D - Function that returns local contribution of the
%                      1D stiffness matrix
%
% INPUT: 
%   xcoord - physical coordinates 
%   diff   - diffusion coefficient  
%   nquad  - number of quadrature points
%
% OUTPUT: 
%   k - local contribution of the stiffness matrix
%

%Get coordinates
x1 = xcoord(1);
x2 = xcoord(2);

%Size of the element
h = x2-x1;

%Quadrature points and weights
[z, w] = lobpts(nquad);

%Derivative of shape functions
dN1 = @(x) -1/h;
dN2 = @(x) 1/h;

%Jacobian
J = h/2;

%Initialize local matrix
k = zeros(2,2);

%Loop through quadrature points 
for i = 1:nquad
    
    %Get the evaluation point
    l = x1+(h/2)*(1+z(i));
        
    %Accumulate the local stiffness matrix   
    k = k+[dN1(l)*dN1(l) dN1(l)*dN2(l);...
         dN2(l)*dN1(l) dN2(l)*dN2(l)]*diff*w(i)*abs(J);
end

end

