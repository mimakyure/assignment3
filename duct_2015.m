
while  U_change > U_change_max
  n_count=n_count+1;

  for i=2:nhx-1
    for j=2:nhy-1
      Unew(i,j)= U(i,j)-dt*(P(i+1,j)-P(i-1,j))/(2*hx)...
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
  U_change=max(max(abs((Unew-U)/dt)));
  U_change;
  U=Unew;

  for i=2:nhx-1
    for j=2:nhy-1
      Vnew(i,j)= V(i,j)-dt*(P(i,j+1)-P(i,j-1))/(2*hy)...
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
  V=Vnew;

  for i=2:nhx-1
    for j=2:nhy-1
      div(i,j)=(U(i+1,j)-U(i-1,j))/(2*hx)...
      +(V(i,j+1)-V(i,j-1))/(2*hy);
    end
  end

  div_sum=sum(sum(abs(div)))/((nhx-2)*(nhy-2));

  resid_pc=1;
  Pc(:,:)=0;
  p_count=0;

  while resid_pc > resid_pc_max & p_count < 100;
    p_count=p_count+1;
    for i=2:nhx-1
      for j=2:nhy-1
        residual(i,j)= -1/(hx*hx)*(Pc(i+1,j)-2.*Pc(i,j)+Pc(i-1,j))...
                                 -1/(hy*hy)*(Pc(i,j+1)-2.*Pc(i,j)+Pc(i,j-1))...
                               +1/dt*div(i,j);
        Pc(i,j)= (1/(-2/(hx*hx)-2/(hy*hy))*residual(i,j))*relx_pc+Pc(i,j);
      end
    end

    Pc(1,:)=Pc(2,:);
    Pc(:,1)=Pc(:,2);
    Pc(:,nhy)=Pc(:,nhy-1);
    Pc(nhx,:)=0;

    resid_pc=sum(sum(abs(residual)))/((nhx-2)*(nhy-2));

  end

  P=P+Pc;

  for i=2:nhx-1
    for j=2:nhy-1
      U(i,j)= U(i,j)-dt*(Pc(i+1,j)-Pc(i-1,j))/(2*hx);
      V(i,j)= V(i,j)-dt*(Pc(i,j+1)-Pc(i,j-1))/(2*hy);
    end
  end

  U(nhx,:)=U(nhx-1,:);
  U(1,2:nhy-1)=-U(2,2:nhy-1)+2;
  U(2:nhx,1)=-U(2:nhx,2);
  U(2:nhx,nhy)=-U(2:nhx,nhy-1);

  V(nhx,:)=V(nhx-1,:);
  V(1,:)=-V(2,:);
  V(:,1)=-V(:,2);
  V(:,nhy)=-V(:,nhy-1);
  U_max=max(max(U));

  if round(n_count/10)*10 == n_count 
    n_count
    div_sum
    U_change
    U_max
  end

end

U_change
