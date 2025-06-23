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
  S_hat = 100;
  %%
   c_I = [230/256 139/256 2/256];
  c_I_face = [239/256 197/256 127/256];
  c_S = [4/256 113/256 181/256];
  c_S_face = [129/256 184/256 218/256];
  c_N = [0 0 0];
  %%%%%
  %%
  b=[exp(1) exp(1.1) exp(1.2)];
  n=3;
  parms=parm_table(n,beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);
  %%
  fig=figure
  fig.Position = [400 400 250 130]
  std=sqrt(variance);
  x = (-5 * std(1):0.01:5 * std(1)) + load(1);  %// Plotting range
  y = exp(- 0.5 * ((x - load(1)) / std(1)) .^ 2) / (std(1) * sqrt(2 * pi));
  plot(x, y, 'LineWidth', 1,'Color',[14/256 204/256 163/256]);hold on
  area(x, y, 'FaceColor', [14/256 204/256 163/256], 'EdgeColor', 'none','FaceAlpha',0.3); % Fill with a blue shade
  x = (-5 * std(2):0.01:5 * std(2)) + load(2);  %// Plotting range
  y = exp(- 0.5 * ((x - load(2)) / std(2)) .^ 2) / (std(2) * sqrt(2 * pi));
  plot(x, y, 'LineWidth', 1,'Color',[11/256 158/256 191/256]);hold on
  area(x, y, 'FaceColor', [11/256 158/256 191/256], 'EdgeColor', 'none','FaceAlpha',0.2); % Fill with a blue shade
  x = (-5 * std(3):0.01:5 * std(3)) + load(3);  %// Plotting range
  y = exp(- 0.5 * ((x - load(3)) / std(3)) .^ 2) / (std(3) * sqrt(2 * pi));
  plot(x, y, 'LineWidth', 1,'Color',[20/256 93/256 227/256]);hold on
  area(x, y, 'FaceColor', [20/256 93/256 227/256], 'EdgeColor', 'none','FaceAlpha',0.2); % Fill with a blue shade
  xlim([-3 15]);ylim([0 0.25])
  xlabel('Log pathogen load'); xlim([-3 12])
    ylabel('Density');
    %lgd=legend([l1 l2 l3],["$b=e$","$b=e^{1.1}$","$b=e^{1.2}$"],"Interpreter","latex",'Location','northwest')
    %set(lgd, 'Box', 'off', 'Color', 'none','ItemTokenSize', [10, 8]);
    grid on;
    set(gca,"Fontsize",14,'FontName', 'Times')
    exportgraphics(fig,'figure/b_dist.pdf','Resolution',600)
%%
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
  S_hat = 100;
  %%
  r_max=[0.1 0.3 0.99];
  n=3;
  parms=parm_table(n,beta,r_max,a,b,mu_0,sigma_0,K,lambda,d_z,l,d_0,r,gamma);
  [S,I,P,Q,Z,load,variance,R0,stability]=Bd_system(S_hat,parms);

%%
fig=figure
  fig.Position = [400 400 250 130]
  std=sqrt(variance);
  x = (-5 * std(1):0.01:5 * std(1)) + load(1);  %// Plotting range
  y = exp(- 0.5 * ((x - load(1)) / std(1)) .^ 2) / (std(1) * sqrt(2 * pi));
  plot(x, y, 'LineWidth', 1,'Color',[14/256 204/256 163/256]);hold on
  area(x, y, 'FaceColor', [14/256 204/256 163/256], 'EdgeColor', 'none','FaceAlpha',0.3); % Fill with a blue shade
  x = (-5 * std(2):0.01:5 * std(2)) + load(2);  %// Plotting range
  y = exp(- 0.5 * ((x - load(2)) / std(2)) .^ 2) / (std(2) * sqrt(2 * pi));
  plot(x, y, 'LineWidth', 1,'Color',[11/256 158/256 191/256]);hold on
  area(x, y, 'FaceColor', [11/256 158/256 191/256], 'EdgeColor', 'none','FaceAlpha',0.2); % Fill with a blue shade
  x = (-5 * std(3):0.01:5 * std(3)) + load(3);  %// Plotting range
  y = exp(- 0.5 * ((x - load(3)) / std(3)) .^ 2) / (std(3) * sqrt(2 * pi));
  plot(x, y, 'LineWidth', 1,'Color',[20/256 93/256 227/256]);hold on
  area(x, y, 'FaceColor', [20/256 93/256 227/256], 'EdgeColor', 'none','FaceAlpha',0.2); % Fill with a blue shade
  xlim([-3 15]);ylim([0 0.25])
  xlabel('Log pathogen load'); xlim([-3 12])
    ylabel('Density');
    %lgd=legend([l1 l2 l3],["$b=e$","$b=e^{1.1}$","$b=e^{1.2}$"],"Interpreter","latex",'Location','northwest')
    %set(lgd, 'Box', 'off', 'Color', 'none','ItemTokenSize', [10, 8]);
    grid on;
    set(gca,"Fontsize",14,'FontName', 'Times')
    exportgraphics(fig,'figure/rmax_dist.pdf','Resolution',600)













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
