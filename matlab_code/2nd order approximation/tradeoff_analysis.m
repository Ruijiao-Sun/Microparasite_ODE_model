clc;clear;
%% parameters
   beta = 1.8e-4;
    d_0 = 0.001; %background mortality rate
  r_max = 0.664;%
   mu_0 = 3.382;%log(0) %oscillation when mu_0 = 5
sigma_0 = 2.708;
      K = exp(12);
 lambda = 0.052;%
    d_z = 4.658;%
      l = 0.003207;%
      a = 2e-5;%1.5e-5
      b = exp(1.06);
      r = 0.013;%daily maximum recruitment rate
  gamma = 0.033; %density dependent strength
      n = 100;
  S_hat = -log(d_0 / r) / gamma;
  %%
  %%
min_value = 0.01;
max_value = 1;
rmax_range = linspace(min_value,max_value,n);
%%
fig=figure
fig.Position = [400 400 200 200]
b = exp(1.06);
  parms=parm_table(n,beta,rmax_range,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
 h1= plot(rmax_range,R0,"Color",[0.8 0.8 0.8],'LineWidth',1.5);hold on
 plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.8 0.8 0.8],'LineWidth',1.5)
 b = exp(1.1);
  parms=parm_table(n,beta,rmax_range,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
  h2=plot(rmax_range,R0,"Color",[0.5 0.5 0.5],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.5 0.5 0.5],'LineWidth',1.5);
b = exp(1.14);
  parms=parm_table(n,beta,rmax_range,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
h3=plot(rmax_range,R0,"k",'LineWidth',1.5);hold on
plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"k--",'LineWidth',1.5);
ylabel ('R_0');xlabel("r_{max}");xlim([min_value,max_value]); %ylim([2 8])
box off
%legend([h1, h2, h3], "$b=e^{1.06}$","$b=e^{1.1}$","$b=e^{1.14}$",'interpreter','latex')
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/tradeoff_b.pdf','Resolution',600)

%%
%%
fig=figure
fig.Position = [400 400 200 200]
d_z=5.5;
  parms=parm_table(n,beta,rmax_range,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
 h1= plot(rmax_range,R0,"Color",[0.8 0.8 0.8],'LineWidth',1.5);hold on
 plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.8 0.8 0.8],'LineWidth',1.5)
 d_z=6;
  parms=parm_table(n,beta,rmax_range,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
  h2=plot(rmax_range,R0,"Color",[0.5 0.5 0.5],'LineWidth',1.5);hold on
  plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"--","Color",[0.5 0.5 0.5],'LineWidth',1.5);
d_z=6.2;
  parms=parm_table(n,beta,rmax_range,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
h3=plot(rmax_range,R0,"k",'LineWidth',1.5);hold on
plot(rmax_range(R0==max(R0))*ones(1,100),linspace(0,max(R0),100),"k--",'LineWidth',1.5);
ylabel ('R_0');xlabel("r_{max}");xlim([min_value,max_value]); ylim([0.5 2.5]);xlim([0 1])
box off
%legend([h1, h2, h3], "$b=e^{1.06}$","$b=e^{1.1}$","$b=e^{1.14}$",'interpreter','latex')
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/tradeoff_dz.pdf','Resolution',600)


%%
















%%
function parms=parm_table(N,beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma)
parms = [beta(:)   .* ones(N, 1), ...
    r_max(:).* ones(N, 1), ...
    a(:)       .* ones(N, 1), ...
    b(:)       .* ones(N, 1), ...
    mu_0(:)    .* ones(N, 1), ...
    sigma_0(:) .* ones(N, 1), ...
    K(:)       .* ones(N, 1), ...
    lambda(:)  .* ones(N, 1), ...
    d_z(:)     .* ones(N, 1), ...
    l(:)       .* ones(N, 1), ...
    d_0(:)     .* ones(N, 1), ...
    r(:)       .* ones(N, 1), ...
    gamma(:)   .* ones(N, 1)];
end
