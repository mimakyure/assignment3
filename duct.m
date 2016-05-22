
% AMME5202
% Semester 1, 2016
% Matthew Imakyure

if exist('OCTAVE_VERSION', 'builtin') ~= 0;
  page_screen_output(0);
  page_output_immediately(1);
end

fprintf('v2\n');

%%
% code timing -----------------------------------------------------------------
% [iter, total time, v/u-calc time, p-calc time, other-calc time]
%
% v1, Gauss-Seidel, No optimisations
% [8795, 5032.4, 37.434, 1686.3, 12.974]
% [8795, 6797.9, 37.472, 1715.5, 12.444]
%
% v2, Gauss-Jacobi, Vectorized operations
% [9347, 7204.7, 10.556, 391.47, 4.7131]
% [9347, 7608.0, 10.418, 388.13, 4.6793]
%
% v2, Gauss-Jacobi, Pressure matrix solver
% 


%%
% given parameters ------------------------------------------------------------
nu     = 0.001;
rho    = 1.0;
Uin    = 1.0;
len    = 4.0;
height = 0.1;


%%
% solver settings -------------------------------------------------------------

% discretisation controls
% hx and hy to align with domain size
dt = 1e-3;
hx = 0.02;
hy = 0.002;

% stopping criteria
U_change_max = 1e-3;
resid_pc_max = 1e-3;

% relaxation factor (not used)
relax = 1;

% might be useful
Cr = 1.5*Uin*dt/hx;
VN = 1*dt/hx^2;
fprintf('Cr = %1.2g\nVN = %1.2g\n', Cr, VN);
fprintf('Cr + 2VN = %1.2g\n', Cr + 2*VN);
fprintf('4VN(1 - VN) - Cr^2 = %1.2g\n', 4*VN*(1 - VN)- Cr^2);

%%
% initialize variables --------------------------------------------------------

% set number of mesh nodes in x and y directions
% ghost cells at edges to maintain boundary conditions
nhx = len/hx + 2;
nhy = height/hy + 2;

% for plotting results
midy = round(nhy/2);
xn = -hx/2:hx:len+hx/2;
yn = -hy/2:hy:height+hy/2;

runs = 3;
timerTimes = zeros(runs, 5);
for timerI = 1:runs

  % number of time step iterations needed to find solution
  n_count = 0;

  % mesh matrixes to compute U velocity, V velocity, and pressure
  U    = zeros(nhx,nhy);
  Unew = zeros(nhx,nhy);
  V    = zeros(nhx,nhy);
  Vnew = zeros(nhx,nhy);
  P    = zeros(nhx,nhy);
  div  = zeros(nhx,nhy);
  residual = zeros(nhx,nhy);

  % indexes for mesh internal nodes
  i = 2:nhx-1;
  j = 2:nhy-1;

  % set initial change in velocity to get in loop
  U_change = 10;

  %%
  % calculate the solution ------------------------------------------------------
  tTime = tic;
  vTot = 0;
  pTot = 0;
  oTot = 0;
  while  U_change > U_change_max
    vTime = tic;
    n_count = n_count + 1;

    %%
    % solve for velocity

    % calculate intermediate U velocity
    % second order central x & y, first order upwind x, first order central
    % central y
    Unew(i,j) = U(i,j) + ...
      + nu*dt*(1/(hx*hx)*(U(i+1,j) - 2*U(i,j) + U(i-1,j)) ...
             + 1/(hy*hy)*(U(i,j+1) - 2*U(i,j) + U(i,j-1))) ...
      - dt*U(i,j)/(hx).*(U(i,j) - U(i-1,j)) ...
      - dt*V(i,j)/(2*hy).*(U(i,j+1) - U(i,j-1));

    % boundary conditions
    Unew(nhx,:)       =  Unew(nhx-1,:);       % zero gradient at outlet
    Unew(1,2:nhy - 1) = -Unew(2,2:nhy-1) + 2; % inlet velocity of 1
    Unew(2:nhx,1)     = -Unew(2:nhx,2);       % zero velocity at wall
    Unew(2:nhx,nhy)   = -Unew(2:nhx,nhy-1);   % zero velocity at wall

    % calculate intermediate V velocity (similar to U equations)
    Vnew(i,j) = V(i,j)...
      + nu*dt*(1/(hx*hx)*(V(i+1,j) - 2*V(i,j) + V(i-1,j)) ...
             + 1/(hy*hy)*(V(i,j+1) - 2*V(i,j) + V(i,j-1))) ...
      - dt*U(i,j)/(hx).*(V(i,j) - V(i-1,j)) ...
      - dt*V(i,j)/(2*hy).*(V(i,j+1) - V(i,j-1));

    % boundary conditions
    Vnew(nhx,:) =  Vnew(nhx-1,:); % zero gradient at outlet
    Vnew(1,:)   = -Vnew(2,:);     % zero velocity at inlet
    Vnew(:,1)   = -Vnew(:,2);     % zero velocity at wall
    Vnew(:,nhy) = -Vnew(:,nhy-1); % zero velocity at wall
    vTot = vTot + toc(vTime);

    %%
    % solve Poisson equation for pressure using matrix solver
    % compute divergence at each node
    pTime = tic;
    resid_pc = 1;
    p_count = 0;

    % compute divergence at each node using first order central scheme
    div(i,j) = (Unew(i+1,j) - Unew(i-1,j))/(2*hx) ... 
      + (Vnew(i,j+1) - Vnew(i,j-1))/(2*hy);

    % calculate average divergence over all nodes (for tracking)
    div_sum = sum(sum(abs(div)))/((nhx - 2)*(nhy - 2));

    % solve Poisson equation for pressure as heat equation
    while resid_pc > resid_pc_max && p_count < 100;
      p_count = p_count + 1;

      % error
      residual(i,j) = (rho/dt)*div(i,j) ...
        - 1/(hx*hx)*(P(i+1,j) - 2.*P(i,j) + P(i-1,j)) ...
        - 1/(hy*hy)*(P(i,j+1) - 2.*P(i,j) + P(i,j-1));

      P(i,j) = (1/(-2/(hx*hx) - 2/(hy*hy))*residual(i,j))*relax + P(i,j);

      % boundary conditions
      % zero pressure gradient at walls
      P(1,:)   = P(2,:);
      P(:,1)   = P(:,2);
      P(:,nhy) = P(:,nhy-1);
      P(nhx,:) = 0;

      resid_pc = sum(sum(abs(residual)))/((nhx - 2)*(nhy - 2));

    end
    pTot = pTot + toc(pTime);

    %%
    % correct velocity based on pressure results
    % central approximation for pressure
    oTime = tic;
    Unew(i,j) = Unew(i,j) - (dt/rho)*(P(i+1,j) - P(i-1,j))/(2*hx);
    Vnew(i,j) = Vnew(i,j) - (dt/rho)*(P(i,j+1) - P(i,j-1))/(2*hy);

    % boundary conditions as above
    Unew(nhx,:)     = Unew(nhx-1,:);
    Unew(1,2:nhy-1) = -Unew(2,2:nhy-1) + 2;
    Unew(2:nhx,1)   = -Unew(2:nhx,2);
    Unew(2:nhx,nhy) = -Unew(2:nhx,nhy-1);

    Vnew(nhx,:) = Vnew(nhx-1,:);
    Vnew(1,:)   = -Vnew(2,:);
    Vnew(:,1)   = -Vnew(:,2);
    Vnew(:,nhy) = -Vnew(:,nhy-1);

    U_change = max(max(abs((Unew - U)/dt)));
    U_max = max(max(Unew(i,:)));

    % pointer swap
    U = Unew;
    V = Vnew;

    % print a status update
    if round(n_count/10)*10 == n_count
      fprintf('%d) div=%1.2e, U_change=%1.2e, U_max=%1.2e\n', ...
        n_count, div_sum, U_change, U_max);
    end
    oTot = oTot + toc(oTime);

  end
  tTot = tTot + toc(tTime);
  fprintf('%d) div=%1.2e, U_change=%1.2e, U_max=%1.2e\n', ...
    n_count, div_sum, U_change, U_max);

  fprintf('Total elapsed time is: %1.6f seconds.\n', tTot);
  fprintf('Velocity calculation elapsed time is: %1.6f seconds (%0.1f%%).\n', ...
    vTot, vTot/tTot*100);
  fprintf('Pressure calculation elapsed time is: %1.6f seconds (%0.1f%%).\n', ...
    pTot, pTot/tTot*100);
  fprintf('Other calculation elapsed time is: %1.6f seconds (%0.1f%%).\n', ...
    oTot, oTot/tTot*100);

  timerTimes(timerI,:) = [n_count, tTot, vTot, pTot, oTot];
end

%%
% plots

% development of max velocity
figure(1);
plot(xn, U(:,midy));
xlim([0 4]);
title('Velocity at Duct Centreline');
xlabel('X-Position (m)');
ylabel('U Velocity (m/s)');

% velocity profiles
figure(2);
profiles = bsxfun(@minus, U(:,j), min(U(:,j), [], 2));
x_profiles = bsxfun(@plus, xn', profiles);
plot(x_profiles(floor(size(x_profiles, 1)/100)*[1,2,4,8,16,32,64,100],:), ...
  yn(j));
title('Velocity Profile Development');
xlabel('X-Position (m)');
ylabel('Y-Position (m)');

