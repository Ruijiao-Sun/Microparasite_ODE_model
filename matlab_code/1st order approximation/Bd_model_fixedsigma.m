function dxdt = Bd_model_fixedsigma(t,x,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma)
         S = x(1); %susceptible individuals
         I = x(2); %infected $individuals
         P = x(3); %log pathogen load
         Z = x(4); %free living zoospores at natural scale
         N = S + I;
         dxdt = zeros(4,1);   %creat empty vector to store variables
         %calculate mean pathogen load in infected individuals
          if I<=0
          %update bd dynamics
             dxdt(1) = Recruit(N,r,gamma)- d_0*S - beta*Z*S; %S
             dxdt(2) = beta*Z*S - l*I; %I
             dxdt(3) = beta*Z*S*mu_0;
             dxdt(4) = -d_z*Z;%Z
          else
             mu = P / I;
             alpha_t = mu +  log(b)*sigma_F^2/2;
             theta_t = mu +  sigma_F^2/2; %derivative quantity

             %update bd dynamics
             dxdt(1) = Recruit(N,r,gamma)- d_0*S - beta*Z*S + l*I; %S
             dxdt(2) = beta*Z*S - l*I - load_death(alpha_t,d_0,a,b)*I; %I
             dxdt(3) = beta*Z*S*mu_0 ...
                     + I*Growth(theta_t,r_max,K) ...
                     - I*(mu*(l+load_death(alpha_t,d_0,a,b)) + sigma_F^2*a*log(b)*b^alpha_t);
             dxdt(4) = -d_z*Z + lambda*exp(theta_t)*I;%Z
          end

end
