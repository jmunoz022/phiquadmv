function dir = BC(nx,ny)
%
% BC - Function that returns the Dirichlet nodes 
%                   
% INPUT: 
%   nx - number of nodes in the x direction
%   ny - number of nodes in the x direction
%
% OUTPUT: 
%   dir - vector containing the Dirichlet nodes
%

dir = zeros(1,2*ny+2*nx-4);
dir(1:nx) = 1:nx;
for i = 1:ny-2
   dir(nx+2*i-1) = i*nx+1;
   dir(nx+2*i) = (i+1)*nx;
end
dir(nx+2*(ny-2)+1:end) = nx*(ny-1)+1:nx*ny;

end