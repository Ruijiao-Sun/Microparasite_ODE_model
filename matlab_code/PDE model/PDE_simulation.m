clc;clear;
%%
results3=main_script()
%%
% Solver
function results = solve_infection_dynamics(tmax, dt, xmin, xmax, nx, ...
    params, S0, Z0)
    % Setup meshgrid
    x = linspace(xmin, xmax, nx);
    dx = x(2) - x(1);
    nt = floor(tmax/dt);
    t = linspace(0, tmax, nt);
    
    % 
    S = zeros(1, nt);
    Z = zeros(1, nt);
    I = zeros(nx, nt);
    
    % Initial conditions with normal distribution
    S(1) = S0;
    Z(1) = Z0;
    I(:,1) = normpdf(x, params.mu_0, params.sigma_0)';
    
    % Normalize initial distribution to have specified total infected
    % population
    I(:,1) = params.I0_total * I(:,1) / trapz(x, I(:,1));
    
    % 
    p0_x = normpdf(x, params.mu_0, params.sigma_0)';  % Using same distribution for p0
    d_x = load_death(x, params.d_0, params.a, params.b)';
    
    % Time step
    for n = 1:nt-1
        % advection 
        G_x = Growth(x, params.r_max, params.K)';
        I_temp = I(:,n);
        
        % difference 
        for j = 2:nx-1
            if G_x(j) > 0
                flux = (G_x(j)*I(j,n) - G_x(j-1)*I(j-1,n))/dx;
            else
                flux = (G_x(j+1)*I(j+1,n) - G_x(j)*I(j,n))/dx;
            end
            I_temp(j) = I_temp(j) - dt * flux;
        end
        
        % Integrals
        recovery_integral = trapz(x, params.l * I(:,n));
        shedding_integral = trapz(x, params.lambda * exp(x)' .* I(:,n));
        
        % Update S
        N = S(n) + trapz(x, I(:,n));
        S(n+1) = S(n) + dt * (...
            Recruit(N, params.r, params.gamma) - ...
            params.d_0*S(n) - params.beta*Z(n)*S(n) + recovery_integral...
        );
        
        % 4. Update Z
        Z(n+1) = Z(n) + dt * (...
            -params.d_z*Z(n) + shedding_integral...
        );
        
        % 5. Update I
        I(:,n+1) = I_temp + dt * (...
            params.beta*Z(n)*S(n)*p0_x - (params.l + d_x).*I_temp...
        );
        
        % Ensure non-negativity
        S(n+1) = max(0, S(n+1));
        Z(n+1) = max(0, Z(n+1));
        I(:,n+1) = max(0, I(:,n+1));
    end
    
    % Wrap results
    results.t = t;
    results.x = x;
    results.S = S;
    results.I = I;
    results.Z = Z;
end

% Figure
function results=main_script()
    % Define parameters
    params.tmax = 10000;
    params.dt = 0.1;
    params.xmin = -10;
    params.xmax = 12;
    params.nx = 200;
    
    % Model parameters
    params.d_0 = 0.001;     % baseline death rate
    params.a = 2e-5;      % death load parameter
    params.b = exp(1.06);       % death load parameter
    params.beta = 1.8e-4;   % transmission rate
    params.l = 0.003207;       % recovery rate
    params.d_z = 4.658;     % zoospore loss rate
    params.lambda = 0.052;  % shedding rate
    params.r = 0.013;       % recruitment rate
    params.gamma = 0.033; % density dependence
    params.r_max = 0.95;   % maximum growth rate
    params.K = exp(12);      % carrying capacity
    params.mu_0 = 3.382;      % Initial load mean
 params.sigma_0 = 2.708;   % Initial load std
    
    % Initial distribution parameters
    params.I0_total = 1; % total initial infected population
    
    % Initial conditions
    S0 = 10;
    Z0 = 10;
    
    % Run simulation
    results = solve_infection_dynamics(params.tmax, params.dt, ...
        params.xmin, params.xmax, params.nx, ...
        params, S0, Z0);
    
    % Visualize results
    visualize_results(results, params);
end

function visualize_results(results, params)  
    % Plot S(t)
    fig=figure
    fig.Position = [400 400 200 250]
    plot(results.t, results.S, 'LineWidth', 2);
    xlabel('Time');
    ylabel('Susceptible Population');
    title('S(t)');
    grid on;
    set(gca,"Fontsize",14,'FontName', 'Times')
    
    % Plot Z(t)
    fig=figure
    fig.Position = [400 400 200 250]
    total_infected = trapz(results.x, results.I, 1);
    total_pop = results.S + total_infected;
    plot(results.t, total_pop, 'LineWidth', 2);
    xlabel('Time');
    ylabel('Total Population');
    title('Total Population vs Time');
    grid on;
    set(gca,"Fontsize",14,'FontName', 'Times')
    
    % Plot I(x,t) as heatmap
    fig=figure
    fig.Position = [400 400 400 250]
    imagesc(results.t, results.x, results.I);
    xlabel('Time');
    ylabel('Log pathogen load');
    %title('I(x,t)');
    c=colorbar;
    axis xy;  % Put zero at bottom
    c.Label.String = 'Density';
    colormap(cmocean('thermal'));
    set(gca,"Fontsize",14,'FontName', 'Times')

    % Plot final distribution
    fig=figure
    fig.Position = [400 400 300 150]
    plot(results.x, results.I(:,end), 'LineWidth', 1.5,'Color',[230/256 139/256 2/256]); hold on
    area(results.x, results.I(:,end), 'FaceColor', [239/256 197/256 127/256], 'EdgeColor', 'none','FaceAlpha',0.5); % Fill with a blue shade
    hold on;
    mean_value = sum(results.x'.* results.I(:,end),"all") / sum(results.I(:,end))
    xline(mean_value, '--r', 'LineWidth', 0.1, 'DisplayName', 'Mean');% Dashed red line for mean
    xlabel('Log pathogen load');
    ylabel('Density');
    grid on;
    set(gca,"Fontsize",14,'FontName', 'Times')
end