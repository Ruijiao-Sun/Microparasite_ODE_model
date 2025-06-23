function equ = Bd_equlibrium(beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma)
        
       syms x y I  
       % x = mu_t; P/I
       % y = simga_t^2; Q/I-(P/I)^2

       alpha = x + log(b)*y / 2;
       theta = x + y / 2;

        eqn1 = mu_0 * (l + d_0 + a * b^alpha) ...
             + r_max * (1 - exp(theta) / K) ...
             - x * (l+d_0 + a * b^alpha) ...
             - y * a * log(b)*b^alpha == 0;

        eqn2 = (mu_0^2+sigma_0^2) * (l + d_0 + a * b^alpha)...
             + 2*((x + y)*r_max * (1 - exp(theta) / K)-y*r_max)...
             - (x^2 + y) * (l + d_0 + a * b^alpha) ...
             - 2 * log(b) * y * alpha * a * b^alpha == 0;

        [mu_equ, sigma_equ2] = vpasolve([eqn1, eqn2], [x,y]);

        S_equ = (l+d_0+a*b^(mu_equ + log(b)*sigma_equ2/2))*d_z/(beta*lambda*exp(mu_equ + sigma_equ2/2));

        eqn3 = r*(S_equ+I)*exp(-gamma*(S_equ+I))-d_0*S_equ - beta*S_equ*lambda*I*exp(mu_equ + sigma_equ2/2)/d_z + l*I ==0;
        I_equ = vpasolve(eqn3,I);
        P_equ = I_equ*mu_equ;
        Q_equ = (sigma_equ2 + mu_equ^2) * I_equ;
        Z_equ = lambda*I_equ*exp(mu_equ + sigma_equ2/2)/d_z;

        equ=double([S_equ;I_equ;P_equ;Q_equ;Z_equ]);
end