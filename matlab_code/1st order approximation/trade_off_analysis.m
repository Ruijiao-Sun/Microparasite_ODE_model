clc;clear;
%%
%equalibrium as a function of r_max, maximum within-host pathogen growth rate
beta = 1.8e-4;
d_0 = 0.001; %background mortality rate
r_max = 0.664;%
sigma_F = 2.708;% 2.4 --    value 2.3 is stable oscillation
mu_0 = 3.382;%log(0) %oscillation when mu_0 = 5
sigma_0 = 2.708;
K = exp(12);
lambda = 0.052;%
d_z = 4.658;%
l = 0.003207;%
a = 1e-5;%1.5e-5
b = exp(1.3);
r = 0.013;%daily maximum recruitment rate
gamma = 0.033; %density dependent strength
n=100;
%% solve for S_hat, disease free equilibrium
S_hat = -log(d_0 / r) / gamma;
%%
min_value = 0.01;
max_value = 2;
rmax_range = linspace(min_value,max_value,n);
%%
fig=figure
fig.Position = [400 400 200 200]
b = exp(1.06);
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
 h1= plot(rmax_range,R0,"Color",[0.8 0.8 0.8],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.8 0.8 0.8],'LineWidth',1.5)
b = exp(1.1);
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
  h2=plot(rmax_range,R0,"Color",[0.5 0.5 0.5],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.5 0.5 0.5],'LineWidth',1.5);
b = exp(1.14);
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
h3=plot(rmax_range,R0,"k",'LineWidth',1.5);hold on
plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"k--",'LineWidth',1.5);
ylabel ('R_0');xlabel("r_{max}");xlim([min_value,max_value]); ylim([2 8])
box off
%legend([h1, h2, h3], "$b=e^{1.06}$","$b=e^{1.1}$","$b=e^{1.14}$",'interpreter','latex')
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/tradeoff_b.pdf','Resolution',600)


%%
fig=figure;
fig.Position = [400 400 200 200];
b=exp(1.14);
sigma_F = 2.5;
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
  h1=plot(rmax_range,R0,"Color",[0.8 0.8 0.8],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.8 0.8 0.8],'LineWidth',1.5)
sigma_F = 2.8;
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
  h2=plot(rmax_range,R0,"Color",[0.5 0.5 0.5],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.5 0.5 0.5],'LineWidth',1.5);
sigma_F = 3;
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
h3=plot(rmax_range,R0,"k",'LineWidth',1.5);hold on
plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"k--",'LineWidth',1.5);
ylabel ('R_0');xlabel("r_{max}");xlim([min_value,max_value]); ylim([1 3.5])
box off
%legend([h1, h2, h3], '\sigma_F=2.5', '\sigma_F=2.8', '\sigma_F=3');
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/tradeoff_sigma.pdf','Resolution',600)

%%
%we only have infection from environmental zoospore pool, no direct
%infecion from infected individuals, equilibria mean pathogen load does not
%depend on dz or lambda, in R0 calculation, d_z does not interact with
%equilibria mean load, so for a given d_z value, R0 is always maximized at
%the same r_max.   Day 2002, ecology letters
%no direct or death-mediated transmission Bonhoeffer, S., Lenski, R.E. & Ebert, D. (1996). The curse of the pharaoh: the evolution of virulence in pathogens with long living propagules. Proc. Royal Soc. London,


b=exp(1.06);
sigma_F = 2.5;
fig=figure;
fig.Position = [400 400 200 200];
d_z=5.5;
%lambda=0.01
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
 h1= plot(rmax_range,R0,"Color",[0.8 0.8 0.8],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.8 0.8 0.8],'LineWidth',1.5)

d_z=6;
%lambda=0.05
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
 h2= plot(rmax_range,R0,"Color",[0.5 0.5 0.5],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.5 0.5 0.5],'LineWidth',1.5);
  rmax_range(R0==max(R0))*ones(1,100)
d_z=6.2;
%lambda=0.1
  parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Z,load,R0]=Bd_system(S_hat,parms);
h3=plot(rmax_range,R0,"k",'LineWidth',1.5);hold on
plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"k--",'LineWidth',1.5);
rmax_range(R0==max(R0))*ones(1,100)
ylabel ('R_0');xlabel("r_{max}");xlim([min_value,max_value]);ylim([4,7]);
box off
% legend([h1, h2, h3], "d_z=5.5","d_z=6","d_z=6.2")
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/tradeoff_dz.pdf','Resolution',600)





















%%

function parms=parm_table(N,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma)
parms = [beta    * ones(N, 1), ...
    r_max(:).* ones(N, 1), ...
    sigma_F * ones(N, 1), ...
    a       * ones(N, 1), ...
    b       * ones(N, 1), ...
    mu_0    * ones(N, 1), ...
    K       * ones(N, 1), ...
    lambda  * ones(N, 1), ...
    d_z     * ones(N, 1), ...
    l       * ones(N, 1), ...
    d_0     * ones(N, 1), ...
    r       * ones(N, 1), ...
    gamma   * ones(N, 1)];

end

