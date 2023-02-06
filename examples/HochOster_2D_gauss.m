%Hochburck-Osterman Equation in 2D+time with phiquadmv and Gauss quadrature
%with Exponential RK methods up to order 3
mkdir(pwd,'Results_HochOster_gauss')
mypath = fullfile(pwd,'Results_HochOster_gauss');

tol = 1.0e-14;

rmax = 5; %Maximum elements in each space direction: 2^r
mmax = 6; %Maximum time step size: 2^m

%Data of the problem
eps = 1;
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;
T = 1;
uexact = @(x,y,t) exp(t).*x.*(1-x).*y.*(1-y);
u0 = @(x,y) x.*(1-x).*y.*(1-y);
source = @(x,y,t) uexact(x,y,t)+2*exp(t).*(x.*(1-x)+y.*(1-y))-1./(1+exp(2*t).*(x.*(1-x).*y.*(1-y)).^2);

for r = 5:rmax 

    %Initialize vectors for the results
    Errors1 = zeros(mmax,1);
    Errors2 = zeros(mmax,1);
    Errors3 = zeros(mmax,1);
    times_Eulerkron = zeros(mmax,1);
    times_RK2kron = zeros(mmax,1);
    times_RK3kron = zeros(mmax,1);
    taum=zeros(mmax,1);
    nsteps=zeros(mmax,1);
    
    %Meshes, Parameters and BC
    nelx = 2^r;
    nely = 2^r;
    xsol = linspace(x1,x2,nelx+1);
    ysol = linspace(y1,y2,nely+1);
    [nx,ny,nel,nnode,coord,nodes] = parameters(xsol,ysol);
    dirx = [1 nx];
    diry = [1 ny];
    nquad = 2; %quadrature points for the matrices (Lobatto quadrature)
    nquadf = 2; %quadrature points for the source term (Gauss quadrature)
    
    %1D matrices
    [Ax,Mx] = FEM_matrices_1D(eps,dirx,xsol,nelx,nquad);
    nnx = size(Ax,1);
    M = kron(Mx,Mx); 
    normA = 2*norm(Ax,'inf');
    
    %Initial condition vector
    dir = BC(nx,ny);
    [xmat,ymat] = meshgrid(xsol,ysol);
    U0 = u0(xmat,ymat);
    U0 = reshape(U0',[],1);
    U0(dir) = [];
    dimx = size(U0,1);
    
    UexactT = uexact(xmat,ymat,T);
    UexactT = reshape(UexactT',[],1);
    UexactT(dir) = [];

    for m = 1:mmax      
        %Time grid and step size
        steps = 2^m;
        nsteps(m) = steps;
        tau = T/steps;
        taum(m) = tau;
        t = 0:tau:T;
        
        %Assembly the linear source term
        fprintf('Assembly of the RHS... \n\n')
        flin = zeros(dimx,steps+1);
        flin12 = zeros(dimx,steps);
        flin13 = zeros(dimx,steps);
        flin23 = zeros(dimx,steps);
        for i = 1:steps
            flin(:,i) = source_term(source,t(i),M,nel,nnode,coord,nodes,nquadf,dir);
            flin12(:,i) = source_term(source,t(i)+tau/2,M,nel,nnode,coord,nodes,nquadf,dir);
            flin13(:,i) = source_term(source,t(i)+tau/3,M,nel,nnode,coord,nodes,nquadf,dir);
            flin23(:,i) = source_term(source,t(i)+2*tau/3,M,nel,nnode,coord,nodes,nquadf,dir);
        end
        flin(:,end) = source_term(source,t(end),M,nel,nnode,coord,nodes,nquadf,dir);

        %%%Exponential Euler method%%%
        fprintf('Running Exponential Euler for r=%d and %d time steps\n\n',r,steps)
        Ukron = zeros(dimx,steps+1);
        Ukron(:,1) = U0;
        tic
        [npts,~,sc,~] = setup_quadrature(tau*normA, tol, 1, 'gauss');
        for i = 1:steps
            f1 = 1./(1+Ukron(:,i).^2)+flin(:,i);
            Ukron(:,i+1) = Ukron(:,i)+tau*phiquadmv_q(1,-tau*Ax,f1-actionA(Ax,Ukron(:,i),nnx),tau*normA,sc,npts); 
        end
        times_Eulerkron(m) = toc;
        Errors1(m) = norm(UexactT-Ukron(:,end),'inf');
    
        %RK2 method%%%
        fprintf('Running RK2 for r=%d and %d time steps\n\n',r,steps)
        Ukron = zeros(dimx,steps+1);
        Ukron(:,1) = U0;
        tic
        [npts,~,sc,~] = setup_quadrature(tau*normA, tol, 1, 'gauss');
        [npts12,~,sc12,~] = setup_quadrature(tau*normA/2, tol, 1, 'gauss');
        for i=1:steps
            %stage 1
            f1kron = 1./(1+Ukron(:,i).^2)+flin(:,i);
            %stage 2
            Ukrons1 = Ukron(:,i)+(tau/2)*phiquadmv_q(1,-tau*Ax/2,f1kron-actionA(Ax,Ukron(:,i),nnx),tau*normA/2,sc12,npts12); 
            f2kron = 1./(1+Ukrons1.^2)+flin12(:,i);
            %compute solution
            Ukron(:,i+1) = Ukron(:,i)+tau*phiquadmv_q(1,-tau*Ax,f2kron-actionA(Ax,Ukron(:,i),nnx),tau*normA,sc,npts); 
        end
        times_RK2kron(m) = toc;
        Errors2(m) = norm(UexactT-Ukron(:,end),'inf');
            
        %third-order method%%%
        fprintf('Running RK3 for r=%d and %d time steps\n\n',r,steps)
        Ukron = zeros(dimx,steps+1);
        Ukron(:,1) = U0;
        tic
        [npts,~,sc,~] = setup_quadrature(tau*normA, tol, 1, 'gauss');
        [npts13,~,sc13,~] = setup_quadrature(tau*normA/3, tol, 1, 'gauss');
        [npts23_1,~,sc23_1,~] = setup_quadrature(2*tau*normA/3, tol, 1, 'gauss');
        [npts23_2,~,sc23_2,~] = setup_quadrature(2*tau*normA/3, tol, 2, 'gauss');
        [npts2,~,sc2,~] = setup_quadrature(tau*normA, tol, 2, 'gauss');
        for i = 1:steps
           %stage 1
           f1kron = 1./(1+Ukron(:,i).^2)+flin(:,i);
           %stage 2
           Ukrons1 = Ukron(:,i)+(tau/3)*phiquadmv_q(1,-tau*Ax/3,f1kron-actionA(Ax,Ukron(:,i),nnx),tau*normA/3,sc13,npts13); 
           f2kron = 1./(1+Ukrons1.^2)+flin13(:,i);
           %stage 3
           Ukrons2 = Ukron(:,i)+(4*tau/3)*phiquadmv_q(2,-2*tau*Ax/3,-f1kron+f2kron,2*tau*normA/3,sc23_2,npts23_2)+(2*tau/3)*phiquadmv_q(1,-2*tau*Ax/3,f1kron-actionA(Ax,Ukron(:,i),nnx),2*tau*normA/3,sc23_1,npts23_1);
           f3kron = 1./(1+Ukrons2.^2)+flin23(:,i);
           %compute solution
           Ukron(:,i+1) = Ukron(:,i)+tau*phiquadmv_q(1,-tau*Ax,f1kron-actionA(Ax,Ukron(:,i),nnx),tau*normA,sc,npts)+tau*phiquadmv_q(2,-tau*Ax,(3/2)*(-f1kron+f3kron),tau*normA,sc2,npts2);
        end
        times_RK3kron(m) = toc;       
        Errors3(m) = norm(UexactT-Ukron(:,end),'inf');
    end

    %Save errors and computational times
    Tab_times = table(nsteps,times_Eulerkron,times_RK2kron,times_RK3kron,'VariableNames',{'Nsteps','Times_EulerKron','Times_RK2Kron','Times_RK3Kron'});
    writetable(Tab_times,fullfile(mypath,['TimesHochOster_Kron_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
    
    Tab_errors = table(taum,Errors1,Errors2,Errors3,'VariableNames',{'tau','Err_EulerKron','Err_RK2Kron','Err_RK3Kron'});
    writetable(Tab_errors,fullfile(mypath,['ErrHochOster_Kron_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
end

%Plot error convergence
figure
loglog(taum,Errors1,taum,Errors2,taum,Errors3)
legend('Euler','RK2','RK3','Location','SouthEast')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)


%Function that returns a sigle action
function y = phiquadmv_q(q,A,b,normA,sc,npts)
    Y = phiquadmv(q,A,b,normA,sc,'gauss',npts);
    y = Y(:,q); 
end

%Function to compute the action Ab
function y = actionA(Ax,b,nnx)
    y = reshape(Ax*reshape(b,nnx,nnx),nnx*nnx,1)+reshape(reshape(b,nnx,nnx)*Ax',nnx*nnx,1);
end