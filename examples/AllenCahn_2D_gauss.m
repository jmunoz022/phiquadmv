%Allen-Cahn Equation in 2D+time with phiquadmv and Gauss quadrature
%with Exponential Euler method

tol = 1.0e-14;

%Data of the problem
eps = 1;
N = 7;
alpha = 0.75;
theta = @(x,y) atan((y-0.5)./(x-0.5)).*(x>0.5)+(pi+atan((y-0.5)./(x-0.5))).*(x<=0.5);
u0 = @(x,y) tanh((0.25+0.1*cos(N*theta(x,y))-sqrt((x-0.5).^2+(y-0.5).^2))/(sqrt(2)*alpha));
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;
T = 2e-2;
delta = 0.01; %coefficient of the non-linear term 

%Number of time steps, time step size and time 
steps = 200;
tau = T/steps;
t = 0:tau:T;

% Meshes, Parameters and BC
nquad = 2;
r = 9;
nelx = 2^r+1;
nely = 2^r+1;
xsol = linspace(x1,x2,nelx+1);
ysol = linspace(y1,y2,nely+1);
[nx,ny,nel,nnode,coord,nodes] = parameters(xsol,ysol);
dirx = [];
diry = [];

%1D matrices
Ax = FEM_matrices_1D(eps,dirx,xsol,nelx,nquad);
nnx = size(Ax,1);
normA = 2*norm(Ax,'inf');

%Initial condition vector
[xmat,ymat] = meshgrid(xsol,ysol);
U0 = u0(xmat,ymat);
U0 = reshape(U0',[],1);
dimx = size(U0,1);

%%%Exponential Euler method%%%
fprintf('Running Exponential Euler for r=%d and %d time steps\n\n',r,steps)
Ukron = zeros(dimx,steps+1);
Ukron(:,1) = U0;
tic
[npts,~,sc,~] = setup_quadrature(tau*normA, tol, 1, 'gauss');
for i = 1:steps
    f1 = (Ukron(:,i)-Ukron(:,i).^3)/delta^2;
    Ukron(:,i+1) = Ukron(:,i)+tau*phiquadmv_q(1,-tau*Ax,f1-actionA(Ax,Ukron(:,i),nnx),tau*normA,sc,npts); 
end
times_Eulerkron = toc

%Plot solution 
figure
xlabel({'$x$'},'interpreter','latex')
ylabel({'$y$'},'interpreter','latex')
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)
for i = 1:steps+1
    Umat = reshape(Ukron(:,i),[nelx+1,nely+1])';
    h = pcolor(xsol,ysol,Umat);
    set(h, 'EdgeColor', 'none')
    colorbar, set(colorbar,'TickLabelInterpreter','latex')
    str = sprintf('Time = %.5f s',t(i));
    set(title(str),'Interpreter','Latex');
    axis([0 1 0 1 -1 1])
    pause(0.01)
end


%Function that returns a sigle action
function y = phiquadmv_q(q,A,b,normA,sc,npts)
    Y = phiquadmv(q,A,b,normA,sc,'gauss',npts);
    y = Y(:,q); 
end

%Function to compute the action Ab
function y = actionA(Ax,b,nnx)
    y = reshape(Ax*reshape(b,nnx,nnx),nnx*nnx,1)+reshape(reshape(b,nnx,nnx)*Ax',nnx*nnx,1);
end
    






