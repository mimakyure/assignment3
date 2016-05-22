
% AMME5202
% Semester 1, 2016
% Matthew Imakyure

if exist('OCTAVE_VERSION', 'builtin') ~= 0;
  page_screen_output(0);
  page_output_immediately(1);
end

fprintf('v3\n');


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
bicg_max     = 1e-3;
bicg_iter    = 100;

% relaxation factor (not used)
relax = 1;

% might be useful
Cr = 1.5*Uin*dt/hx;
VN = 1*dt/hx^2;
fprintf('Cr = %1.2g\n', Cr);
fprintf('Cr + 2VN = %1.2g\n', Cr + 2*VN); % < 1
fprintf('VN = %1.2g\n', VN); % < 1/2
fprintf('4VN(1 - VN) - Cr^2 = %1.2g\n', 4*VN*(1 - VN)- Cr^2); % > 0


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

  % mesh matrices to compute U velocity, V velocity, and pressure
  U    = zeros(nhx,nhy);
  Unew = zeros(nhx,nhy);
  V    = zeros(nhx,nhy);
  Vnew = zeros(nhx,nhy);
  P    = zeros(nhx,nhy);
  Pnew = zeros((nhx-2)*(nhy-2),1);
  div  = zeros(nhx,nhy);
  residual = zeros(nhx,nhy);

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

  % P(1,:) = P(2,:) -> P^{x-} = P; zero gradient at inlet
  % P gets P^{x-} coefficient
  A(1) = A(1) + hx^-2;

  % P(nhx,:) = 0 -> P^{x+} = 0; zero pressure at outlet ghost cell
  % P^{x+} gone at outlet

  % expand into system of equations to apply remainder of conditions
  A = kron(eye(Aj), A); % ghost cells remain zero doing this

  % add coefficients for y+/- pressures at +/-nhx, matrix becomes pentadiagonal
  A = A + spdiags(ones(Ai*Aj,1)*[hy^-2 hy^-2], [-Ai Ai], Ai*Aj, Ai*Aj);

  % P(:,1)   = P(:,2)     -> P^{y-} = P; zero gradient at lower wall
  % P(:,nhy) = P(:,nhy-1) -> P^{y+} = P; zero gradient at upper wall
  E = sparse([1:Ai,Ai*(Aj-1)+1:Ai*Aj], [1:Ai,Ai*(Aj-1)+1:Ai*Aj], ...
    ones(2*Ai,1)*hy^-2, Ai*Aj, Ai*Aj);
  A = A + E; % P gets P^{y+} coefficient at top and P^{y-} coefficent at bottom


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
    vTot = vTot + toc(vTime);


    %%
    % solve Poisson equation for pressure using matrix solver
    pTime = tic;
    resid_pc = 1;
    p_count = 0;

    % compute divergence at each node using central difference
    div(i,j) = (Unew(i+1,j) - Unew(i-1,j))/(2*hx) ... 
      + (Vnew(i,j+1) - Vnew(i,j-1))/(2*hy);

    % calculate average divergence over all nodes (for tracking)
    div_sum = sum(sum(abs(div)))/((nhx - 2)*(nhy - 2));

    % calculate source term trimming off ghost cells
    f = rho/dt*div(2:end-1,2:end-1);

    % use iterative matrix solver to resolve pressure
    [Pnew, flag] = bicgstab(A, f(:), bicg_max, bicg_iter);
    P = sparse(nhx, nhy);
    P(2:end-1,2:end-1) = reshape(Pnew, nhx-2, nhy-2); 
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

