clc;clear;
% %% Parameters
%    beta = 1.8e-4;
%     d_0 = 0.001; %background mortality rate
%   r_max = 0.664;%
% sigma_F = 2.708;% 2.4 --    value 2.3 is stable oscillation
%    mu_0 = 3.382;%log(0) %oscillation when mu_0 = 5
% sigma_0 = 2.708;
%       K = exp(12);
%  lambda = 0.052;%
%     d_z = 4.658;%
%       l = 0.003207;%
%       a = 2e-5;%1.5e-5;
%       b = exp(1.06);
%       r = 0.013;%daily maximum recruitment rate
%   gamma = 0.033; %density dependent strength
%       n = 100;
% S_hat = -log(d_0 / r) / gamma;

%% Parameters
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
      a = 2e-5;%1.5e-5;
      b = exp(1.1);
      r = 0.013;%daily maximum recruitment rate
  gamma = 0.033; %density dependent strength
      n = 100;
S_hat = -log(d_0 / r) / gamma;
  %% colors for plot
  c_I = [230/256 139/256 2/256];
  c_I_face = [239/256 197/256 127/256];
  c_S = [4/256 113/256 181/256];
  c_S_face = [129/256 184/256 218/256];
  c_N = [0 0 0];
%%
%equalibrium as a function of sigma
min_value = 1;
max_value = 5;
sigma_F_range = linspace(min_value,max_value,n);
parms=parm_table(n,beta,r_max,sigma_F_range,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
%% figures
fig=figure
fig.Position = [400 400 200 200]
plot(sigma_F_range(stability==1),load(stability==1),"k",'LineWidth',1.5);hold on
plot(sigma_F_range(stability==0),load(stability==0),"-.k",'LineWidth',1.5);hold on
x_stable=sigma_F_range(stability==1);y_stable=load(stability==1);
scatter(x_stable(1:10:end),y_stable(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('Mean log pathogen load');xlabel("$\sigma_F$","Interpreter","latex");xlim([min_value,max_value])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/sigma_load.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
N_stable=S(stability==1)+I(stability==1);S_stable=S(stability==1);I_stable=I(stability==1);
plot(sigma_F_range(stability==1),N_stable,'k','LineWidth',1.5);hold on
scatter(x_stable(1:10:end),N_stable(1:10:end),60,'MarkerEdgeColor',c_N,'MarkerFaceColor',c_N,'LineWidth',.5)
plot(sigma_F_range(stability==1),S(stability==1),'LineWidth',1.5,"Color",c_S);
scatter(x_stable(1:10:end),S_stable(1:10:end),60,'MarkerEdgeColor',c_S,'MarkerFaceColor',c_S_face,'LineWidth',.5)
plot(sigma_F_range(stability==1),I(stability==1),'LineWidth',1.5,"Color",c_I);
scatter(x_stable(1:10:end),I_stable(1:10:end),60,'MarkerEdgeColor',c_I,'MarkerFaceColor',c_I_face,'LineWidth',.5)
plot(sigma_F_range(stability==0),S(stability==0)+I(stability==0),'-.k','LineWidth',1.5);hold on
plot(sigma_F_range(stability==0),S(stability==0),'-.','LineWidth',1.5,"Color",c_S);
plot(sigma_F_range(stability==0),I(stability==0),"-.",'LineWidth',1.5,"Color",c_I);
ylabel ('Population abundance');xlabel("$\sigma_F$","Interpreter","latex");xlim([min_value,max_value])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%exportgraphics(fig,'figure/sigma_pop.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
x_stable=sigma_F_range(stability==1);
plot(sigma_F_range(stability==1),I_stable./N_stable,'k','LineWidth',1.5);hold on
scatter(x_stable(1:10:end),I_stable(1:10:end)./N_stable(1:10:end),60,'MarkerEdgeColor',c_N,'MarkerFaceColor',c_N,'LineWidth',.5)
ylabel ('Prevalence');xlabel("$\sigma_F$","Interpreter","latex");xlim([1,5]);%ylim([0.05 0.08])
box off;ylim([0 0.06])
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/sigma_prevalence.pdf','Resolution',600)
%%
%% explore effects of b on prevalence as a function of sigma
fig=figure
fig.Position = [400 400 200 200]
b_temp=exp(1.12)
parms=parm_table(n,beta,r_max,sigma_F_range,a,b_temp,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
plot(sigma_F_range,I./(S+I),"Color",[15/256 16/256 18/256],'LineWidth',2);hold on
b_temp=exp(1.1)
parms=parm_table(n,beta,r_max,sigma_F_range,a,b_temp,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
plot(sigma_F_range,I./(S+I),"Color",[114/256,130/256,147/256],'LineWidth',2);
b_temp=exp(1.08)
parms=parm_table(n,beta,r_max,sigma_F_range,a,b_temp,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
plot(sigma_F_range,I./(S+I),"Color",[152/256,180/256,183/256],'LineWidth',2);hold on
ylabel ('Prevalence');xlabel("$\sigma_F$","Interpreter","latex");xlim([1,5]);%ylim([0.05 0.08])
box off;ylim([0 0.06]);xlim([1 4])
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%%
fig=figure
fig.Position = [400 400 200 200]
R0_stable=R0(stability==1);
plot(sigma_F_range(stability==1),R0(stability==1),"k",'LineWidth',1.5);hold on
plot(sigma_F_range(stability==0),R0(stability==0),"-.k",'LineWidth',1.5);
scatter(x_stable(1:10:end),R0_stable(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('$R_0$',"Interpreter","latex");xlabel("$\sigma_F$","Interpreter","latex");xlim([min_value,max_value]);%ylim([23.75 23.85])
box off
set(gca,"tickdir",'out',"Fontsize",14,"Fontname","times")
exportgraphics(fig,'figure/sigma_R0.pdf','Resolution',600)
%%
%equalibrium as a function of a, load-dependent mortality curve
min_value = 1e-5;
max_value = 2.5e-5;
a_range = linspace(min_value,max_value,n);
parms=parm_table(n,beta,r_max,sigma_F,a_range,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
%% figures
fig=figure
fig.Position = [400 400 200 200]
plot(a_range,load,"k",'LineWidth',1.5);hold on
scatter(a_range(1:10:end),load(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('Mean log pathogen load');xlabel("a");xlim([min_value,max_value])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%%
fig=figure
fig.Position = [400 400 200 200]
plot(a_range,S+I,"k",'LineWidth',1.5);hold on
scatter(a_range(1:10:end),S(1:10:end)+I(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
plot(a_range,S,'LineWidth',1.5,"Color",c_S);
scatter(a_range(1:10:end),S(1:10:end),60,'MarkerEdgeColor',c_S,'MarkerFaceColor',c_S_face,'LineWidth',.5)
plot(a_range,I,'LineWidth',1.5,"Color",c_I);
scatter(a_range(1:10:end),I(1:10:end),60,'MarkerEdgeColor',c_I,'MarkerFaceColor',c_I_face,'LineWidth',.5)
ylabel ('Population abundance');xlabel("a");xlim([min_value,max_value])
%legend("Total","Susceptible","Infected")
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%%
fig=figure
fig.Position = [400 400 200 200]
plot(a_range,R0,"k",'LineWidth',1.5);hold on
plot(a_range,R0,"-.k",'LineWidth',1.5);
scatter(a_range(1:10:end),R0(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('$R_0$',"Interpreter","latex");xlabel("a");xlim([min_value,max_value]);%ylim([23.75 23.85])
box off
set(gca,"tickdir",'out',"Fontsize",14,"Fontname","times")
%%
%equalibrium as a function of b, load-dependent mortality curve
min_value = exp(0.98);
max_value = exp(1.1);
b_range = linspace(min_value,max_value,n);
parms=parm_table(n,beta,r_max,sigma_F,a,b_range,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
%% figures
fig=figure
fig.Position = [400 400 200 200]
plot(b_range,load,"k",'LineWidth',1.5);hold on
scatter(b_range(1:10:end),load(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('Mean log pathogen load');xlabel("b");xlim([min_value,max_value])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/b_load.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
plot(b_range,S+I,"k",'LineWidth',1.5);hold on
scatter(b_range(1:10:end),S(1:10:end)+I(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
plot(b_range,S,'LineWidth',1.5,"Color",c_S);
scatter(b_range(1:10:end),S(1:10:end),60,'MarkerEdgeColor',c_S,'MarkerFaceColor',c_S_face,'LineWidth',.5)
plot(b_range,I,'LineWidth',1.5,"Color",c_I);
scatter(b_range(1:10:end),I(1:10:end),60,'MarkerEdgeColor',c_I,'MarkerFaceColor',c_I_face,'LineWidth',.5)
ylabel ('Population abundance');xlabel("b");xlim([min_value,3])
%legend("Total","Susceptible","Infected")
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%exportgraphics(fig,'figure/b_pop.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
plot(b_range,I./(S+I),"k",'LineWidth',1.5);hold on
scatter(b_range(1:10:end),I(1:10:end)./(S(1:10:end)+I(1:10:end)),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('Prevalence');xlabel("b");xlim([min_value,3])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%exportgraphics(fig,'figure/b_prevalence.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
plot(b_range,R0,"k",'LineWidth',1.5);hold on
plot(b_range,R0,"-.k",'LineWidth',1.5);
scatter(b_range(1:10:end),R0(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('$R_0$',"Interpreter","latex");xlabel("b");xlim([min_value,max_value]);%ylim([23.75 23.85])
box off
set(gca,"tickdir",'out',"Fontsize",14,"Fontname","times")
exportgraphics(fig,'figure/b_R0.pdf','Resolution',600)
%%
%equalibrium as a function of r_max, maximum within-host pathogen growth rate
min_value = 0.1;
max_value = 2;
%b=exp(0.97);
rmax_range = linspace(min_value,max_value,n);
parms=parm_table(n,beta,rmax_range,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma);
[S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms);
%% figures
fig=figure
fig.Position = [400 400 200 200]
plot(rmax_range,load,"k",'LineWidth',1.5);hold on
scatter(rmax_range(1:10:end),load(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('Mean log pathogen load');xlabel("$r_{max}$","Interpreter","latex");xlim([min_value,max_value])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
exportgraphics(fig,'figure/rmax_load.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
plot(rmax_range,S+I,"k",'LineWidth',1.5);hold on
scatter(rmax_range(1:10:end),S(1:10:end)+I(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
plot(rmax_range,S,'LineWidth',1.5,"Color",c_S);
scatter(rmax_range(1:10:end),S(1:10:end),60,'MarkerEdgeColor',c_S,'MarkerFaceColor',c_S_face,'LineWidth',.5)
plot(rmax_range,I,'LineWidth',1.5,"Color",c_I);
scatter(rmax_range(1:10:end),I(1:10:end),60,'MarkerEdgeColor',c_I,'MarkerFaceColor',c_I_face,'LineWidth',.5)
ylabel ('Population abundance');xlabel("r_{max}");xlim([min_value,max_value])
%legend("Total","Susceptible","Infected")
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%exportgraphics(fig,'figure/rmax_pop.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
plot(rmax_range,I./(S+I),"k",'LineWidth',1.5);hold on
scatter(rmax_range(1:10:end),I(1:10:end)./(S(1:10:end)+I(1:10:end)),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('Prevalence');xlabel("r_{max}");xlim([min_value,max_value]);ylim([0 0.25])
box off
set(gca,"tickdir",'out',"Fontsize",14,'FontName', 'Times')
%exportgraphics(fig,'figure/rmax_prevalence.pdf','Resolution',600)
%%
fig=figure
fig.Position = [400 400 200 200]
plot(rmax_range,R0,"k",'LineWidth',1.5);hold on
plot(rmax_range,R0,"-.k",'LineWidth',1.5);
scatter(rmax_range(1:10:end),R0(1:10:end),60,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'LineWidth',.5)
ylabel ('$R_0$',"Interpreter","latex");xlabel("r_{max}");xlim([min_value,max_value]);%ylim([23.75 23.85])
box off
set(gca,"tickdir",'out',"Fontsize",14,"Fontname","times")
exportgraphics(fig,'figure/rmax_R0.pdf','Resolution',600)
%%






















%%
function parms=parm_table(N,beta,r_max,sigma_F,a,b,mu_0,K,lambda,d_z,l,d_0,r,gamma)
parms = [beta(:)   .* ones(N, 1), ...
    r_max(:).* ones(N, 1), ...
    sigma_F(:) .* ones(N, 1), ...
    a(:)       .* ones(N, 1), ...
    b(:)       .* ones(N, 1), ...
    mu_0(:)    .* ones(N, 1), ...
    K(:)       .* ones(N, 1), ...
    lambda(:)  .* ones(N, 1), ...
    d_z(:)     .* ones(N, 1), ...
    l(:)       .* ones(N, 1), ...
    d_0(:)     .* ones(N, 1), ...
    r(:)       .* ones(N, 1), ...
    gamma(:)   .* ones(N, 1)];
end
