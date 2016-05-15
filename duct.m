
while  U_change > U_change_max
  n_count=n_count+1;

  for i=2:nhx-1
    for j=2:nhy-1
      Unew(i,j)= U(i,j)...
        +nu*dt*(1/(hx*hx)*(U(i+1,j)-2.*U(i,j)+U(i-1,j))...
               +1/(hy*hy)*(U(i,j+1)-2.*U(i,j)+U(i,j-1)))...
	-dt*U(i,j)/(hx)*(U(i,j)-U(i-1,j))...
        -dt*V(i,j)/(2*hy)*(U(i,j+1)-U(i,j-1));
    end
  end
  Unew(nhx,:)=Unew(nhx-1,:);
  Unew(1,2:nhy-1)=-Unew(2,2:nhy-1)+2;
  Unew(2:nhx,1)=-Unew(2:nhx,2);
  Unew(2:nhx,nhy)=-Unew(2:nhx,nhy-1);

  for i=2:nhx-1
    for j=2:nhy-1
      Vnew(i,j)= V(i,j)...
        +nu*dt*(1/(hx*hx)*(V(i+1,j)-2.*V(i,j)+V(i-1,j))...
               +1/(hy*hy)*(V(i,j+1)-2.*V(i,j)+V(i,j-1)))...
        -dt*U(i,j)/(hx)*(V(i,j)-V(i-1,j))...
        -dt*V(i,j)/(2*hy)*(V(i,j+1)-V(i,j-1));
    end
  end
  Vnew(nhx,:)=Vnew(nhx-1,:);
  Vnew(1,:)=-Vnew(2,:);
  Vnew(:,1)=-Vnew(:,2);
  Vnew(:,nhy)=-Vnew(:,nhy-1);

  for i=2:nhx-1
    for j=2:nhy-1
      div(i,j)=(Unew(i+1,j)-Unew(i-1,j))/(2*hx)...
              +(Vnew(i,j+1)-Vnew(i,j-1))/(2*hy);
    end
  end

  div_sum=sum(sum(abs(div)))/((nhx-2)*(nhy-2));

  resid_pc=1;
  p_count=0;

  while resid_pc > resid_pc_max & p_count < 100;
    p_count=p_count+1;
    for i=2:nhx-1
      for j=2:nhy-1
        residual(i,j)= -1/(hx*hx)*(P(i+1,j)-2.*P(i,j)+P(i-1,j))...
                       -1/(hy*hy)*(P(i,j+1)-2.*P(i,j)+P(i,j-1))...
                       +(rho/dt)*div(i,j);
        P(i,j)= (1/(-2/(hx*hx)-2/(hy*hy))*residual(i,j))*relx_pc+P(i,j);
      end
    end

    P(1,:)=P(2,:);
    P(:,1)=P(:,2);
    P(:,nhy)=P(:,nhy-1);
    P(nhx,:)=0;

    resid_pc=sum(sum(abs(residual)))/((nhx-2)*(nhy-2));

  end

  for i=2:nhx-1
    for j=2:nhy-1
      Unew(i,j)= Unew(i,j)-(dt/rho)*(P(i+1,j)-P(i-1,j))/(2*hx);
      Vnew(i,j)= Vnew(i,j)-(dt/rho)*(P(i,j+1)-P(i,j-1))/(2*hy);
    end
  end

  Unew(nhx,:)=Unew(nhx-1,:);
  Unew(1,2:nhy-1)=-Unew(2,2:nhy-1)+2;
  Unew(2:nhx,1)=-Unew(2:nhx,2);
  Unew(2:nhx,nhy)=-Unew(2:nhx,nhy-1);

  Vnew(nhx,:)=Vnew(nhx-1,:);
  Vnew(1,:)=-Vnew(2,:);
  Vnew(:,1)=-Vnew(:,2);
  Vnew(:,nhy)=-Vnew(:,nhy-1);

  U_change=max(max(abs((Unew-U)/dt)));
  U_max=max(max(Unew));

  % pointer swap
  Utmp = U;
  U = Unew;
  Unew = U;

  Vtmp = V;
  V = Vnew;
  Vnew = V;
  
  if round(n_count/10)*10 == n_count 
    n_count
    div_sum
    U_change
    U_max
  end

end

U_change
