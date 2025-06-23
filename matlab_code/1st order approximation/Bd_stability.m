function stability = Bd_stability(S,I,P,Z,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma)
       Bd_J=Bd_jacobian(S,I,P,Z,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
       stability=all(real(eig(Bd_J)) < 0);
end