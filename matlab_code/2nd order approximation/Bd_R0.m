function R0=Bd_R0(S_hat,I,P,Q,beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0)

F=[0 0 0 beta*S_hat;...
   0 0 0 mu_0*beta*S_hat;...
   0 0 0 (mu_0^2+sigma_0^2)*beta*S_hat;...
   0 0 0 0];

dVIdI = d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2) - I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2);
dVIdP = I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2);
dVIdQ = (a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2)/2;
dVIdZ = 0;

dVPdI = r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1) - I*((P*(d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)))/I^2 - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(Q/I^2 - (2*P^2)/I^3) - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2)*(- P^2/I^2 + Q/I) + (P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2))/I) + (P*(d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)))/I - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(- P^2/I^2 + Q/I) - (I*r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^3 + P/I^2 + Q/(2*I^2)))/K;
dVPdP = I*((d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2))/I + (2*P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b))/I^2 - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(1/I - (P*log(b))/I^2)*(- P^2/I^2 + Q/I) + (P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2))/I) - (I*r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(P/I^2 - 1/I))/K;
dVPdQ = (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I)))/(2*K) - I*((a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b))/I + (a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^3*(- P^2/I^2 + Q/I))/(2*I) - (P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2)/(2*I^2));
dVPdZ = 0;

dVQdI = 2*r_max*(- P^2/I^2 + Q/I) - 2*I*(r_max*(Q/I^2 - (2*P^2)/I^3) + r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1)*(P/I^2 + Q/I^2 - (2*P^2)/I^3) + (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^2 + P/I + Q/I)*(- P^2/I^3 + P/I^2 + Q/(2*I^2)))/K) + 2*r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1)*(- P^2/I^2 + P/I + Q/I) - Q*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2) + 2*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(- P^2/I^2 + Q/I) - 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(Q/I^2 - (2*P^2)/I^3) - 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2)*(- P^2/I^2 + Q/I) - 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2)*(- P^2/I^2 + Q/I);
dVQdP = Q*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2) - 2*I*(r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1)*((2*P)/I^2 - 1/I) + (2*P*r_max)/I^2 + (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(P/I^2 - 1/I)*(- P^2/I^2 + P/I + Q/I))/K) + 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2)*(- P^2/I^2 + Q/I) - (4*P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2))/I + 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(1/I - (P*log(b))/I^2)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(- P^2/I^2 + Q/I);
dVQdQ = d_0 + l + 2*I*(r_max/I + (r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1))/I + (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^2 + P/I + Q/I))/(2*I*K)) + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2) + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(- P^2/I^2 + Q/I) + 2*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2) + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^3*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(- P^2/I^2 + Q/I) + (Q*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2)/(2*I);
dVQdZ = 0;

dVZdI = I*lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^3 + P/I^2 + Q/(2*I^2)) - lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I));
dVZdP = I*lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(P/I^2 - 1/I);
dVZdQ = -(lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I)))/2;
dVZdZ = d_z;



V=[dVIdI dVIdP dVIdQ dVIdZ;...
   dVPdI dVPdP dVPdQ dVPdZ;...
   dVQdI dVQdP dVQdQ dVQdZ;...
   dVZdI dVZdP dVZdQ dVZdZ];

R0=double(max(eig(F/V)));

end
