function [nx,ny,nel,nnode,coord,nodes] = parameters(xsol,ysol)
%
% PARAMETERS - Function that returns different parameters related to 
%              the space discretization 
%
% INPUT: 
%   xsol - spatial grid in the x direction         
%   ysol - spatial grid in the y direction
%
% OUTPUT: 
%   nx    - number of nodes in the x direction 
%   ny    - number of nodes in the y direction 
%   nel   - total number of elements 
%   nnode - total number of nodes 
%   coord - vector containing the physical coordinates
%   nodes - vector numbering the nodes
%

%Parameters
nx = size(xsol,2);
ny = size(ysol,2);
nel = (nx-1)*(ny-1);
nnode = nx*ny;

%Coordinates
coord = zeros(nnode,2);
coord(:,1) = repmat(xsol,1,ny);
r = repmat(ysol,nx,1);
coord(:,2) = r(:);

%Nodes 
nodes = zeros(nel,4);
l = 0;
m = 0;
for k = 1:ny-1
    for i = 1:nx-1
        nodes(i+l,1) = i+m+l;
        nodes(i+l,2) = i+1+m+l;
        nodes(i+l,3) = i+nx+1+m+l;
        nodes(i+l,4) = i+nx+m+l;
    end
    l = l+nx-1;
    m = m+1;
end

end