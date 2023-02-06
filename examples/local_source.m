function fn=local_source(source,nquad,xcoord,ycoord,t0)
%%%Function that calculates the local source vector at a fixed time tn%%%

%source=source function
%nquad=number of quadrature points
%xcoord=x coordinates of the element
%ycoord=y coordinates of the element

%%%coordinates of the element%%%
x1 = xcoord(1);
x2 = xcoord(2);
y1 = ycoord(1);
y4 = ycoord(4);

%%%length of the element%%%
a = x2-x1;        
b = y4-y1; 

%%%Quadrature points and weights%%%
[z,w] = legpts(nquad);

%%%basis functions%%% 
N1 = @(x,y) (1-x)*(1-y)/4;
N2 = @(x,y) (1+x)*(1-y)/4;
N3 = @(x,y) (1+x)*(1+y)/4;
N4 = @(x,y) (1-x)*(1+y)/4;

%%%Change of variables%%%
P = @(x,y) x1*(N1(x,y)+N4(x,y))+x2*(N2(x,y)+N3(x,y));
Q = @(x,y) y1*(N1(x,y)+N2(x,y))+y4*(N3(x,y)+N4(x,y));

%%%Jacobian%%%
J = a*b/4;

%%%Integration%%%
fn = zeros(4,1);
I1 = 0;
I2 = 0;
I3 = 0;
I4 = 0;
for i = 1:nquad
    for j = 1:nquad
        I = w(i)*w(j)*source(P(z(i),z(j)),Q(z(i),z(j)),t0)*abs(J);
        I1 = I1+I*N1(z(i),z(j));
        I2 = I2+I*N2(z(i),z(j));
        I3 = I3+I*N3(z(i),z(j));
        I4 = I4+I*N4(z(i),z(j));
    end
end
fn(1) = I1;
fn(2) = I2;
fn(3) = I3;
fn(4) = I4;

%%%Matlab quadrature%%%
% fn1 = zeros(4,1);
% N1 = @(x,y) (x2-x).*(y4-y)/(a*b);
% N2 = @(x,y) (x-x1).*(y4-y)/(a*b);
% N3 = @(x,y) (x-x1).*(y-y1)/(a*b);
% N4 = @(x,y) (x2-x).*(y-y1)/(a*b);
% fn1(1) = quad2d(@(x,y) source(x,y,t0).*N1(x,y),x1,x2,y1,y4,'AbsTol',1e-15);
% fn1(2) = quad2d(@(x,y) source(x,y,t0).*N2(x,y),x1,x2,y1,y4,'AbsTol',1e-15);
% fn1(3) = quad2d(@(x,y) source(x,y,t0).*N3(x,y),x1,x2,y1,y4,'AbsTol',1e-15);
% fn1(4) = quad2d(@(x,y) source(x,y,t0).*N4(x,y),x1,x2,y1,y4,'AbsTol',1e-15);



end