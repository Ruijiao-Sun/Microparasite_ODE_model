function J = Bd_jacobian(S,I,P,Z,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma)
J11=r*exp(-gamma*(I + S)) - Z*beta - d_0 - gamma*r*exp(-gamma*(I + S))*(I + S);
J12=l + r*exp(-gamma*(I + S)) - gamma*r*exp(-gamma*(I + S))*(I + S);
J13=0;
J14=-S*beta;

J21=Z*beta;
J22=(P*a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b))/I - l - a*b^((log(b)*sigma_F^2)/2 + P/I) - d_0;
J23=-a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b);
J24=S*beta;

J31=Z*beta*mu_0;
J32=I*((P*(d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I)))/I^2 + (P^2*a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b))/I^3 + (P*a*b^((log(b)*sigma_F^2)/2 + P/I)*sigma_F^2*log(b)^2)/I^2) - r_max*(exp(sigma_F^2/2 + P/I)/K - 1) - (P*(d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I)))/I - a*b^((log(b)*sigma_F^2)/2 + P/I)*sigma_F^2*log(b) + (P*r_max*exp(sigma_F^2/2 + P/I))/(I*K);
J33= - I*((d_0 + l + a*b^((log(b)*sigma_F^2)/2 + P/I))/I + (P*a*b^((log(b)*sigma_F^2)/2 + P/I)*log(b))/I^2 + (a*b^((log(b)*sigma_F^2)/2 + P/I)*sigma_F^2*log(b)^2)/I) - (r_max*exp(sigma_F^2/2 + P/I))/K;
J34=S*beta*mu_0;

J41=0;
J42=lambda*exp(sigma_F^2/2 + P/I) - (P*lambda*exp(sigma_F^2/2 + P/I))/I;
J43=lambda*exp(sigma_F^2/2 + P/I);
J44=-d_z;

J=[J11 J12 J13 J14;
   J21 J22 J23 J24;
   J31 J32 J33 J34;
   J41 J42 J43 J44];

end