function dxdt = Bd_model_varysigma(t,x,beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma)
         S = x(1); %susceptible individuals
         I = x(2); %infected $individuals
         P = x(3); %log pathogen load
         Q = x(4);
         Z = x(5); %free living zoospores at natural scale
         N = S + I;
         dxdt = zeros(5,1);   %creat empty vector to store variables
         
         %calculate mean pathogen load in infected individuals
          if I<=0
          %update bd dynamics
             dxdt(1) = Recruit(N,r,gamma)- d_0*S - beta*Z*S; %S
             dxdt(2) = beta*Z*S - l*I; %I
             dxdt(3) = beta*Z*S*mu_0;
             dxdt(4) = (sigma_0^2+mu_0^2)*beta*Z*S;
             dxdt(5) = -d_z*Z;%Z
          else
             mu = P / I;
             sigma_t2=Q/I-mu^2;
             alpha_t = mu +  log(b)*sigma_t2/2;
             theta_t = mu +  sigma_t2/2; %derivative quantity
            
             %update bd dynamics
             dxdt(1) = Recruit(N,r,gamma)- d_0*S - beta*Z*S + l*I; %S
             dxdt(2) = beta*Z*S - l*I - load_death(alpha_t,d_0,a,b)*I; %I
             dxdt(3) = beta*Z*S*mu_0 ...
                     + I*Growth(theta_t,r_max,K) ...
                     - I*(mu*(l+load_death(alpha_t,d_0,a,b)) + sigma_t2*a*log(b)*b^alpha_t);
             dxdt(4) = (sigma_0^2+mu_0^2)*beta*Z*S...
                     + 2*I*((mu+sigma_t2)*Growth(theta_t,r_max,K)-sigma_t2*r_max)...
                     - (sigma_t2+mu^2)*I*(l+load_death(alpha_t,d_0,a,b))...
                     - 2*I*log(b)*sigma_t2*alpha_t*a*b^alpha_t;
             dxdt(5) = -d_z*Z + lambda*exp(theta_t)*I;%Z
          end

end
