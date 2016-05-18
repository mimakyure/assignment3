% AMME5202
% Semester 1, 2016
% Matthew Imakyure

if exist('OCTAVE_VERSION', 'builtin') ~= 0;
  page_screen_output(0);
  page_output_immediately(1);
end

%%
% code timing -----------------------------------------------------------------
tic;
%
% v1, no optimisations, 12.5635 s

%%
% given parameters ------------------------------------------------------------
nu     = 0.001;
rho    = 1;
Uin    = 1.0;
len    = 4;
height = 0.1;


%%
% solver settings -------------------------------------------------------------

% discretisation controls
% hx and hy are adjusted later to align with domain size
dt = 1e-3;
hx = 0.0025;
hy = 0.025;

% stopping criteria
U_change_max = 1e-3;
resid_pc_max = 1e-1;

% relaxation factor (not used)
relax = 1;

% might be useful
Cr = 1.5*Uin*dt/hx;
VN = 1*dt/hx^2;

%%
% initialize variables --------------------------------------------------------

% number of time step iterations needed to find solution
n_count = 0;

% set number of mesh nodes in x and y directions
% ghost cells at edges to maintain boundary conditions
nhx = len/hx + 3;
nhy = height/hy + 3;
midy = round(nhy/2); % for plotting centreline velocity

% mesh matrixes to compute U velocity, V velocity, and pressure
U    = zeros(nhx,nhy);
Unew = zeros(nhx,nhy);
V    = zeros(nhx,nhy);
Vnew = zeros(nhx,nhy);
P    = zeros(nhx,nhy);
div  = zeros(nhx,nhy);

% indexes for mesh internal nodes
i = 2:nhx-1;
j = 2:nhy-1;

% set initial change in velocity to get in loop
U_change = 10;

%%
% calculate the solution ------------------------------------------------------

while  U_change > U_change_max
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

  %%
  % solve Poisson equation for pressure using matrix solver
  % compute divergence at each node
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

    P(i,j) = (1/(-2/(hx*hx) - 2/(hy*hy))*residual(i,j))*relx_pc + P(i,j);

    % boundary conditions
    % zero pressure gradient at walls
    P(1,:)   = P(2,:);
    P(:,1)   = P(:,2);
    P(:,nhy) = P(:,nhy-1);
    P(nhx,:) = 0;

    resid_pc = sum(sum(abs(residual)))/((nhx - 2)*(nhy - 2));

  end

  %%
  % correct velocity based on pressure results
  % central approximation for pressure
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

end
toc


fprintf('%d) div=%1.2e, U_change=%1.2e, U_max=%1.2e\n', ...
  n_count, div_sum, U_change, U_max);

%%
% plots

% node locations to plot results
xn = -hx:hx:len+hx;
yn = -hy:hy:height+hy;

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

