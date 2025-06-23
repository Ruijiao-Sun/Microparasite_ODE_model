function J = Bd_jacobian(S,I,P,Q,Z,beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma)
 J11 = r*exp(-gamma*(I + S)) - Z*beta - d_0 - gamma*r*exp(-gamma*(I + S))*(I + S);
 J12 = l + r*exp(-gamma*(I + S)) - gamma*r*exp(-gamma*(I + S))*(I + S);
 J13 = 0;
 J14 = 0;
 J15 = -S*beta;

 J21 = Z*beta;
 J22 = I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2) - l - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2) - d_0;
 J23 = -I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2);
 J24 = -(a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2)/2;
 J25 = S*beta;

 J31 = Z*beta*mu_0;
 J32 = I*((P*(d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)))/I^2 + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(Q/I^2 - (2*P^2)/I^3) + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2)*(- P^2/I^2 + Q/I) + (P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2))/I) - r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1) - (P*(d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)))/I - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(- P^2/I^2 + Q/I) + (I*r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^3 + P/I^2 + Q/(2*I^2)))/K;
 J33 = (I*r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(P/I^2 - 1/I))/K - I*((d_0 + l + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2))/I - (2*P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b))/I^2 + a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(1/I - (P*log(b))/I^2)*(- P^2/I^2 + Q/I) + (P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2))/I);
 J34 = - I*((a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b))/I + (a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^3*(- P^2/I^2 + Q/I))/(2*I) + (P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2)/(2*I^2)) - (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I)))/(2*K);
 J35 = S*beta*mu_0;

 J41 = Z*beta*(mu_0^2 + sigma_0^2);
 J42 = 2*I*(r_max*(Q/I^2 - (2*P^2)/I^3) + r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1)*(P/I^2 + Q/I^2 - (2*P^2)/I^3) + (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^2 + P/I + Q/I)*(- P^2/I^3 + P/I^2 + Q/(2*I^2)))/K) - 2*r_max*(- P^2/I^2 + Q/I) - 2*r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1)*(- P^2/I^2 + P/I + Q/I) + Q*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2) - 2*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(- P^2/I^2 + Q/I) + 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(Q/I^2 - (2*P^2)/I^3) + 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2)*(- P^2/I^2 + Q/I) + 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(P/I^2 + (log(b)*(Q/I^2 - (2*P^2)/I^3))/2)*(- P^2/I^2 + Q/I);
 J43 = 2*I*(r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1)*((2*P)/I^2 - 1/I) + (2*P*r_max)/I^2 + (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(P/I^2 - 1/I)*(- P^2/I^2 + P/I + Q/I))/K) - Q*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2) - 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(1/I - (P*log(b))/I^2)*(- P^2/I^2 + Q/I) + (4*P*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2))/I - 2*I*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(1/I - (P*log(b))/I^2)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(- P^2/I^2 + Q/I);
 J44 = - d_0 - l - 2*I*(r_max/I + (r_max*(exp(- P^2/(2*I^2) + P/I + Q/(2*I))/K - 1))/I + (r_max*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^2 + P/I + Q/I))/(2*I*K)) - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2) - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2*(- P^2/I^2 + Q/I) - 2*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2) - a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^3*(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*(- P^2/I^2 + Q/I) - (Q*a*b^(P/I + (log(b)*(- P^2/I^2 + Q/I))/2)*log(b)^2)/(2*I);
 J45 = S*beta*(mu_0^2 + sigma_0^2);

 J51 = 0;
 J52 = lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I)) - I*lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(- P^2/I^3 + P/I^2 + Q/(2*I^2));
 J53 = -I*lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I))*(P/I^2 - 1/I);
 J54 = (lambda*exp(- P^2/(2*I^2) + P/I + Q/(2*I)))/2;
 J55 = -d_z;

 J=[J11 J12 J13 J14 J15;
    J21 J22 J23 J24 J25;
    J31 J32 J33 J34 J35;
    J41 J42 J43 J44 J45;
    J51 J52 J53 J54 J55];
end