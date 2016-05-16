% First order Euler scheme time discretisation
% Approximate pressure at the new time step and use Poisson equation
page_screen_output(0);
page_output_immediately(1);

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
length = 4;
height = 0.1;


%%
% solver settings -------------------------------------------------------------

% discretisation controls
dt = 1e-3;
hx = 1e-3;
hy = 1e-3;

% set number of mesh nodes in x and y directions
% ghost cells at edges to maintain boundary conditions
nhx = length/hx + 2;
nhy = height/hy + 2;

% stopping criteria
U_change_max = 1e-3;
resid_pc_max = 1e-1;

% pressure solution relaxation factor
relx_pc = 0.1;



%%
% initialize variables --------------------------------------------------------

% number of time step iterations needed to find solution
n_count = 0;

% set initial change in velocity to get in loop
U_change = 10;

% mesh matrixes to compute U velocity, V velocity, and pressure
U    = zeros(nhx,nhy);
Unew = zeros(nhx,nhy);
V    = zeros(nhx,nhy);
Vnew = zeros(nhy,nhx);
P    = zeros(nhx,nhy);

[nhx, nhy] = size(U); % matrix is reversed with rows as x, columns as y

U(1, :) = Uin; % inlet on left

%%
% calculate the solution ------------------------------------------------------

while  U_change > U_change_max
  n_count = n_count + 1;

  % calculate intermediate U velocity
  for i = 2:nhx-1
    for j = 2:nhy-1
      Unew(i,j) = U(i,j) ...
        + nu*dt*(1/(hx*hx)*(U(i+1,j) - 2.*U(i,j) + U(i-1,j)) ...
               + 1/(hy*hy)*(U(i,j+1) - 2.*U(i,j) + U(i,j-1))) ...
	      - dt*U(i,j)/(hx)*(U(i,j) - U(i-1,j)) ...
        - dt*V(i,j)/(2*hy)*(U(i,j+1) - U(i,j-1));
    end
  end

  Unew(nhx,:)       =  Unew(nhx-1,:);
  Unew(1,2:nhy - 1) = -Unew(2,2:nhy-1) + 2;
  Unew(2:nhx,1)     = -Unew(2:nhx,2);
  Unew(2:nhx,nhy)   = -Unew(2:nhx,nhy-1);

  % calculate intermediate V velocity
  for i = 2:nhx-1
    for j = 2:nhy-1
      Vnew(i,j) = V(i,j)...
        + nu*dt*(1/(hx*hx)*(V(i+1,j) - 2.*V(i,j) + V(i-1,j)) ...
               + 1/(hy*hy)*(V(i,j+1) - 2.*V(i,j) + V(i,j-1))) ...
        - dt*U(i,j)/(hx)*(V(i,j) - V(i-1,j)) ...
        - dt*V(i,j)/(2*hy)*(V(i,j+1) - V(i,j-1));
    end
  end

  Vnew(nhx,:) =  Vnew(nhx-1,:);
  Vnew(1,:)   = -Vnew(2,:);
  Vnew(:,1)   = -Vnew(:,2);
  Vnew(:,nhy) = -Vnew(:,nhy-1);

  % compute divergence at each node
  for i = 2:nhx-1
    for j = 2:nhy-1
      div(i,j) = (Unew(i+1,j) - Unew(i-1,j))/(2*hx) ...
        + (Vnew(i,j+1) - Vnew(i,j-1))/(2*hy);
    end
  end

  % calculate average divergence over all nodes
  div_sum = sum(sum(abs(div)))/((nhx - 2)*(nhy - 2));

  % solve Poisson equation for pressure as heat equation
  resid_pc = 1;
  p_count = 0;

  while resid_pc > resid_pc_max & p_count < 100;
    p_count = p_count + 1;
    for i = 2:nhx-1
      for j = 2:nhy-1

        % error
        residual(i,j) = (rho/dt)*div(i,j) ...
          - 1/(hx*hx)*(P(i+1,j) - 2.*P(i,j) + P(i-1,j)) ...
          - 1/(hy*hy)*(P(i,j+1) - 2.*P(i,j) + P(i,j-1));

        P(i,j) = (1/(-2/(hx*hx) - 2/(hy*hy))*residual(i,j))*relx_pc + P(i,j);
      end
    end

    P(1,:)   = P(2,:);
    P(:,1)   = P(:,2);
    P(:,nhy) = P(:,nhy-1);
    P(nhx,:) = 0;

    resid_pc = sum(sum(abs(residual)))/((nhx - 2)*(nhy - 2));

  end

  % correct velocity based on pressure results
  for i = 2:nhx-1
    for j = 2:nhy-1
      Unew(i,j) = Unew(i,j) - (dt/rho)*(P(i+1,j) - P(i-1,j))/(2*hx);
      Vnew(i,j) = Vnew(i,j) - (dt/rho)*(P(i,j+1) - P(i,j-1))/(2*hy);
    end
  end

  % enforce boundary conditions
  Unew(nhx,:)     = Unew(nhx-1,:);
  Unew(1,2:nhy-1) = -Unew(2,2:nhy-1) + 2;
  Unew(2:nhx,1)   = -Unew(2:nhx,2);
  Unew(2:nhx,nhy) = -Unew(2:nhx,nhy-1);

  Vnew(nhx,:) = Vnew(nhx-1,:);
  Vnew(1,:)   = -Vnew(2,:);
  Vnew(:,1)   = -Vnew(:,2);
  Vnew(:,nhy) = -Vnew(:,nhy-1);

  U_change = max(max(abs((Unew - U)/dt)));
  U_max = max(max(Unew));

  % pointer swap
  Utmp = U;
  U = Unew;
  Unew = U;

  Vtmp = V;
  V = Vnew;
  Vnew = V;
  
  % print a status update
  fprintf('%d) div=%1.2e, U_change=%1.2e, U_max=%1.2e\n', ...
    n_count, div_sum, U_change, U_max);

end

U_change

toc
