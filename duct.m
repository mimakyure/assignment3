
% AMME5202
% Semester 1, 2016
% Matthew Imakyure
%
%

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
nu     = 0.001; % kinematic viscosity
rho    = 1.0;   % density
Uin    = 1.0;   % inlet flow velocity
len    = 4.0;   % duct length
height = 0.1;   % duct height


%%
% solver settings -------------------------------------------------------------

% discretisation controls
% duct length and height evenly divisible by hx and hy repectively
dt = 1e-7;
hx = 0.01;
hy = 0.001;

% stopping criteria
U_change_max = 1e-3;
resid_pc_max = 1e-3;
bicg_max     = 1e-3;
bicg_iter    = 100;


% relaxation factor (not used)
relax = 1;

% might be useful for checking stability
Cr = 1.5*Uin*dt/hx;
VN = 1*dt/hx^2;
fprintf('Cr = %1.2g\n', Cr);
fprintf('Cr + 2VN = %1.2g\n', Cr + 2*VN); % < 1
fprintf('VN = %1.2g\n', VN); % < 1/2
fprintf('4VN(1 - VN) - Cr^2 = %1.2g\n', 4*VN*(1 - VN)- Cr^2); % > 0
pause(3);



%%
% initialize variables --------------------------------------------------------

% count number of iterations needed to find solution
n_count = 0;

% set number of mesh nodes in x and y directions
% duct edges are on edge between internal cell and external ghost cell
nhx = len/hx + 2;
nhy = height/hy + 2;

% for plotting results
midy = round(nhy/2);
xn = -hx/2:hx:len+hx/2;
yn = -hy/2:hy:height+hy/2;

% matrixes to store U velocity, V velocity, and pressure on mesh
U    = zeros(nhx,nhy);
Unew = zeros(nhx,nhy);
V    = zeros(nhx,nhy);
Vnew = zeros(nhx,nhy);
P    = zeros(nhx,nhy);
Pnew = zeros((nhx-2)*(nhy-2),1);
div  = zeros(nhx,nhy);

% indexes for mesh internal nodes
i = 2:nhx-1;
j = 2:nhy-1;

% set initial change in velocity to get in loop
U_change = 10;

%%
% assemble pressure coefficients matrix for poisson equation matrix solver
% ghost cells disappear from conditions and embedding done next so not included
Ai = nhx - 2;
Aj = nhy - 2;

% build main tri-diagonal
A = spdiags(ones(Ai, 1)* ...
  [hx^-2 -2*(hx^-2 + hy^-2) hx^-2], ...
  [-1 0 1], Ai, Ai);

% embed boundary conditions by substituting P for P^{x-}, P^{y-}, P^{y+}
% P(1,:) = P(2,:); zero gradient at inlet, P^{x-} gone at inlet
A(1) = A(1) + hx^-2;

% P(nhx,:) = 0; zero pressure at outlet ghost cell; no substitions  here
% P^{x+} gone at outlet

% expand into system of equations to apply remainder of conditions
A = kron(eye(Aj), A);

% add coefficients for y+/- pressures, matrix becomes pentadiagonal
A = A + spdiags(ones(Ai*Aj,1)*[hy^-2 hy^-2], [-Ai Ai], Ai*Aj, Ai*Aj);

% P(:,1) = P(:,2); P(:,nhy) = P(:,nhy-1); zero gradient at walls
% P^{y+} gone at upper wall, P^{y-} gone at lower wall
E = sparse(Ai*Aj, Ai*Aj);
E(1:Ai,1:Ai) = speye(Ai)*hy^-2;
E(end-(Ai-1):end,end-(Ai-1):end) = speye(Ai)*hy^-2;
A = A + E;


%%
% calculate the solution ------------------------------------------------------

while  U_change > U_change_max
  n_count = n_count + 1;

  %%
  % solve for velocity

  % calculate intermediate U velocity
  % 2nd order central x & y, first order upwind x, 1st order central y
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
  resid_pc = 1;
  p_count = 0;

  % compute divergence at each node using central difference
  div(i,j) = (Unew(i+1,j) - Unew(i-1,j))/(2*hx) ... 
    + (Vnew(i,j+1) - Vnew(i,j-1))/(2*hy);

  % calculate source term trimming off ghost cells
  f = rho/dt*div(2:end-1,2:end-1);

  % calculate average divergence over all nodes (for tracking)
  div_sum = sum(sum(abs(div)))/((nhx - 2)*(nhy - 2));

  % use iterative matrix solver to resolve pressure
  [Pnew, flag] = bicgstab(A, f(:), bicg_max, bicg_iter);
  P = sparse(nhx, nhy);
  P(2:end-1,2:end-1) = reshape(Pnew, nhx-2, nhy-2); 


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
    plot(xn, U(:,midy));
    xlim([0 4]);
    drawnow limitrate;
  end

end
toc


fprintf('%d) div=%1.2e, U_change=%1.2e, U_max=%1.2e\n', ...
  n_count, div_sum, U_change, U_max);

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

