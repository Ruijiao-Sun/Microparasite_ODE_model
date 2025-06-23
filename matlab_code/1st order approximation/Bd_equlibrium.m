function equ = Bd_equlibrium(beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma)
         %R0=Bd_R0(S_hat,I,P,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0)
syms x  
syms I
eqn1 = mu_0 * (l + d_0 + a * b^(x + log(b)*sigma_F^2 / 2)) ...
        + r_max * (1 - exp(x + sigma_F^2 / 2) / K) ...
        - x * (l+d_0 + a * b^(x + log(b)*sigma_F^2 / 2)) ...
        - sigma_F^2 * a * log(b)*b^(x + log(b)*sigma_F^2 / 2) == 0;
mu_equ = vpasolve(eqn1, x);
S_equ = (l+d_0+a*b^(mu_equ + log(b)*sigma_F^2/2))*d_z/(beta*lambda*exp(mu_equ + sigma_F^2/2));
eqn2 = r*(S_equ+I)*exp(-gamma*(S_equ+I))-d_0*S_equ - beta*S_equ*lambda*I*exp(mu_equ + sigma_F^2/2)/d_z + l*I ==0;
I_equ = vpasolve(eqn2,I);
P_equ = I_equ*mu_equ;
Z_equ = lambda*I_equ*exp(mu_equ + sigma_F^2/2)/d_z;
equ =double( [S_equ;I_equ;P_equ;Z_equ]);
end