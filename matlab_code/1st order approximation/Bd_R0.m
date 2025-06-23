function R0=Bd_R0(S_hat,I,P,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0)

F=[0 0 beta*S_hat;...
   0 0 mu_0*beta*S_hat;...
   0 0 0];

dVIdI = d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I) - (P*a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b))/I;%diff(VI,I)
dVIdP = a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b);
dVIdZ = 0;

dVPdI =(P*(d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I)))/I - r_max*(exp(sigma_F^2/2 + P/I)/K - 1) - I*((P*(d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I)))/I^2 + (P^2*a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b))/I^3 - (P*a*b^((log(b)*sigma_F^2)/2 + P/I)*sigma_F^2*log(b)^2)/I^2) - a*b^((log(b)*sigma_F^2)/2 + P/I)*sigma_F^2*log(b) + (P*r_max*exp(sigma_F^2/2 + P/I))/(I*K);
dVPdP =I*((d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I))/I + (P*a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b))/I^2 - (a*b^((log(b)*sigma_F^2)/2 + P/I)*sigma_F^2*log(b)^2)/I) - (r_max*exp(sigma_F^2/2 + P/I))/K;
dVPdZ = 0;

dVZdI = (P*lambda*exp(sigma_F^2/2 + P/I))/I - lambda*exp(sigma_F^2/2 + P/I);
dVZdP = -lambda*exp(sigma_F^2/2 + P/I);
dVZdZ = d_z;

V=[dVIdI dVIdP dVIdZ;...
   dVPdI dVPdP dVPdZ;...
   dVZdI dVZdP dVZdZ];

R0=double(max(eig(F/V)));
end