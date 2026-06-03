
clear; clc; close all; 
folder = 'C:\Users\Malin Storm\Documents\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\Building_All\4buildings_DHS_GRID2';
if ~exist(folder,'dir')
     mkdir(folder)
end

% Add all folders inside ./fcns to the MATLAB path.
% This makes helper functions such as getOmega, getPsi, getIdIndex, etc. available.
addpath(genpath('./fcns'));

% Get root path of the BRCM toolbox/project.
BRCMRootPath = getBRCMRootPath();

% Move MATLAB current folder to the BRCM root folder.
cd(BRCMRootPath);

% load building data
building_data = {'building1_data.mat','building2_data.mat','building3_data.mat','building4_data.mat'};


n_buildings = 4; % Number of buildings
n_days = 1; % Number of simulation days
n_sim = n_days*24*4 + 24; % Total simulation steps plus extra prediction buffer
Ts_hrs = 0.25; % Sampling time [h]

Q_DHS_max = 300000; % Maximum DHS power [W]
Q_grid    = 300000; % Maximum grid power [W]
Qmax = 300000; % Common upper power bound [W]

price_tbl = readtable('strompriser_30d_NO1.xlsx'); % Load grid prices
price = price_tbl{:,2}; % Extract price column
price = price(:); % Ensure column vector
price_15min_grid = repelem(price,4); % Convert hourly prices to 15-min resolution
COP = 3.53; % Heat pump coefficient of performance

price_tbl_dhs = readtable('strompriser_30d_NO1_DHS.xlsx'); % Load DHS prices
price_dhs = price_tbl_dhs{:,2}; % Extract price column
price_dhs = price_dhs(:); % Ensure column vector
price_15min_dhs = repelem(price_dhs,4); % Convert hourly prices to 15-min resolution

%% Storage allocation per building
data_ = cell(1,n_buildings);
A_expand = cell(1,n_buildings);
Bu_expand = cell(1,n_buildings);
Bv_expand = cell(1,n_buildings);
Omega = cell(1,n_buildings);
Psi   = cell(1,n_buildings);
Gamma = cell(1,n_buildings);
model = cell(1,n_buildings);
Aineq_ = cell(1,n_buildings);
bineq_ = cell(1,n_buildings);
x_array = cell(1,n_buildings);
u_array = cell(1,n_buildings);
v_array = cell(1,n_buildings);


%% Build MPC model for each building
for i = 1:n_buildings


    data = load(building_data{i}); % Load saved data for building i
    data_{i} = data; % Store building data for later use

    A = data.A; % State matrix
    Bu = data.Bu; % Control input matrix
    Bv = data.Bv; % Disturbance matrix

    nx = data.nx; % Number of building states
    nu = data.nu; % Number of radiator inputs
    nv = data.nv; % Number of disturbances

    Sq_m_zones = data.Sq_m_zones; % Zone areas [m^2]
    V = data.tank_size.V; % Tank volume [m^3]

    Ts = data.Ts_hrs*3600; % Sampling time [s]
    b = Ts/(1000*4186*V); % Tank temperature coefficient

    A_expand{i} = blkdiag(A,1); % Add tank temperature as extra state
    Bu_expand{i} = [Bu, zeros(nx,2);-b*Sq_m_zones', b, b]; % Add DHS/grid inputs and tank dynamics
    Bv_expand{i} = [Bv; zeros(1,nv)]; % Add zero disturbance effect on tank state

    nx_expand = size(A_expand{i},1); % Number of states including tank
    nu_expand = size(Bu_expand{i},2); % Number of inputs including DHS and grid

    n_pred = data.mpc_tuning.n_pred; % Prediction horizon
    nz_zones = size(Sq_m_zones,1); % Number of zones

    ns_ = nz_zones*n_pred; % Total zone-temperature slacks
    nx_expand_ = nx_expand*n_pred; % Total expanded states over horizon
    nu_expand_ = nu_expand*n_pred; % Total expanded inputs over horizon
    %% Prediction matrices
    Omega{i} = getOmega(A_expand{i},n_pred,'lti'); % Maps current state to predicted states
    Psi{i} = getPsi(A_expand{i},Bu_expand{i},n_pred,'lti'); % Maps future inputs to predicted states
    Gamma{i} = kron(eye(n_pred),Bv_expand{i}); % Maps future disturbances to predicted states

    %% Output matrices
    I = eye(nz_zones); % Identity matrix for zone temperatures
    I_nz = zeros(nz_zones,nx_expand); % Output selection matrix
    I_nz(:,1:nz_zones) = I; % Select zone temperature states
    C = I_nz; % Output matrix, y = C*x
    I_nz_ = kron(eye(n_pred), I_nz); % Output matrix over prediction horizon

  
    %% Bounds 
    T_min = 15; T_max = 30; % Zone temperature limits [°C]
    T_min_wt = 35; T_max_wt = 55; % Tank temperature limits [°C]

    lb_x = -inf*ones(nx_expand_,1); % Lower bounds for expanded states
    ub_x = inf*ones(nx_expand_,1); % Upper bounds for expanded states

    for k = 1:n_pred
        tank_indices = (k-1)*nx_expand + nx_expand; % Tank state index at step k
        lb_x(tank_indices) = T_min_wt; % Lower tank temperature bound
        ub_x(tank_indices) = T_max_wt; % Upper tank temperature bound
    end

    u_min = -data.g(1:2:end); % Lower radiator bounds from saved constraint vector
    u_max = data.g(2:2:end); % Upper radiator bounds from saved constraint vector

    u_min_expand = [u_min;0;0]; % Add lower bounds for DHS and grid
    u_max_expand = [u_max;Qmax;Qmax]; % Add upper bounds for DHS and grid

    lb_u = repmat(u_min_expand,n_pred,1); % Lower control bounds over horizon
    ub_u = repmat(u_max_expand,n_pred,1); % Upper control bounds over horizon

    lb_sigma = zeros(ns_,1); % Lower bounds for zone temperature slacks
    ub_sigma = inf(ns_,1); % Upper bounds for zone temperature slacks

    lb_sigma_tank = zeros(n_pred,1); % Lower bounds for tank slacks
    ub_sigma_tank = inf(n_pred,1); % Upper bounds for tank slacks

    lb = [lb_x; lb_u; lb_sigma; lb_sigma_tank]; % Combined lower bounds
    ub = [ub_x; ub_u; ub_sigma; ub_sigma_tank]; % Combined upper bounds

    %% Temperature constraints
    Aineq_temp = [I_nz_, zeros(ns_,nu_expand_), -eye(ns_), zeros(ns_,n_pred);
             -I_nz_, zeros(ns_,nu_expand_), -eye(ns_), zeros(ns_,n_pred)];

    bineq_temp = [T_max*ones(ns_,1);-T_min*ones(ns_,1)];

    %% Tank constraint
    A_tank = zeros(n_pred,length(lb)); % allocate tank constraint matrix
    b_tank = zeros(n_pred,1); % allocate tank constraint RHS

    for k = 1:n_pred
        % Indices in decision vector z = [x; u; sigma; sigma_tank]
        u_indices = nx_expand_ + (k-1)*nu_expand + (1:nu); % Radiator input indices at step k
        tank_indices = (k-1)*nx_expand + nx_expand; % Tank temperature index at step k
        sigma_tank_indices = nx_expand_ + nu_expand_ + ns_ + k; % Tank slack index at step k

        % Tank inequality:
        % Sq_m_zones'*u_rad - (1/b)*T_tank - sigma_tank <= -T_min_wt/b
        A_tank(k,u_indices) = Sq_m_zones'; % zone areas
        A_tank(k,tank_indices) = -1/b; % Tank temperature coefficient
        A_tank(k,sigma_tank_indices) = -1; % Tank slack coefficient
        b_tank(k) = -T_min_wt/b; % RHS of inequality
    end

     
    %% Inequality constraints

    Aineq = [Aineq_temp; A_tank]; % Combine temperature and tank inequality constraints
    bineq = [bineq_temp; b_tank]; % Combine corresponding right-hand sides

    Aineq_{i} = Aineq; % Store inequality matrix for building i
    bineq_{i} = bineq; % Store inequality RHS for building i
    
    %% Cost weights
    Qy = 1000; %1e-6 ; % Weight on zone temperature tracking error

    T_tank_ref = 40; % Tank temperature reference [°C]
    quad_pen = 24;
    Q_tank = 10; % Tank temperature penalty weight

    E_tank = diag([zeros(nx_expand-1,1);1]); % Matrix selecting tank temperature state
    R = data.mpc_tuning.R; % Input penalty matrix from saved building tuning
    
   % Hessian matrix for the quadratic cost function.
    % Blocks: states including tank, control inputs, zone slacks and tank slacks.
    H = blkdiag(kron(eye(n_pred), C'*Qy*C + Q_tank*E_tank),kron(eye(n_pred), R),quad_pen*eye(ns_),1e6*eye(n_pred));

    % Equality constraint for system dynamics.
    % Slack variables have zero columns because they do not affect state evolution.
    A_eq = [eye(nx_expand_),-Psi{i},zeros(nx_expand_,ns_),zeros(nx_expand_,n_pred)];
    
    % Building model
    model{i}.Q = sparse(H);
    model{i}.A = sparse([Aineq_{i}; A_eq]);
    model{i}.sense = [repmat('<',size(Aineq_{i},1),1); repmat('=',size(A_eq,1),1)];
    model{i}.lb = lb;
    model{i}.ub = ub;
    
    % allocating storage
    x_array{i} = zeros(nx_expand,n_sim+1);
    u_array{i} = zeros(nu_expand,n_sim);
   
    % init x-vector 
    x_array{i}(:,1) = [22*ones(nx,1); 45];  % > T_min
end


%% SIMULATION

for kk = 1:(n_sim - n_pred) % simulation loop

    % MPC for each building
    for i = 1:n_buildings % Solve one MPC problem per building

        data = data_{i}; % Get data for building i

        n_pred = data.mpc_tuning.n_pred; % Prediction horizon for building i
        nx_expand = size(A_expand{i},1); % Number of states including tank
        nu_expand = size(Bu_expand{i},2); % Number of inputs including DHS and grid

        nx_expand_ = nx_expand*n_pred; % Total expanded states over horizon
        nu_expand_ = nu_expand*n_pred; % Total expanded inputs over horizon

        %% Predictions
        v_data = data.v_array; % Disturbance data for building i
        weather_idx = data.weather_idx; % Weather disturbance indices

        v_end = min(kk + n_pred - 1, size(v_data,2)); % Last available disturbance index
        v_pred = v_data(:,kk:v_end); % Disturbance prediction over available horizon

        while size(v_pred,2) < n_pred % Pad forecast to full prediction horizon
            v_pred(:,end+1) = v_pred(:,end); % Repeat last available disturbance
        end

        % Forecast noise
        w_temp = 0.2*randn(1,size(v_pred,2)); % Ambient temperature noise
        w_sol  = 15*randn(1,size(v_pred,2)); % Solar radiation noise

        v_pred(weather_idx.Tamb,:) = v_pred(weather_idx.Tamb,:) + w_temp; % Add temperature noise
        v_pred(weather_idx.SolE,:) = max(0, v_pred(weather_idx.SolE,:) + w_sol); % Add east solar noise
        v_pred(weather_idx.SolS,:) = max(0, v_pred(weather_idx.SolS,:) + w_sol); % Add south solar noise

        v_pred = v_pred(:); % Stack disturbance prediction into one vector

        %% Equality constraint
        B_eq = Omega{i}*x_array{i}(:,kk) + Gamma{i}*v_pred; % Dynamic RHS for building i

        nz_zones = size(data.Sq_m_zones,1); % Number of zones in building i
        ns_ = nz_zones*n_pred; % Total zone-temperature slacks over horizon

        n_dec = nx_expand_ + nu_expand_ + ns_ + n_pred; % Total number of decision variables

        model{i}.rhs = [bineq_{i}; B_eq]; % Update constraint RHS
        f = zeros(n_dec,1); % Allocate linear cost vector

        %% Output matrix
        C = [eye(nz_zones), zeros(nz_zones, nx_expand-nz_zones)]; % Select zone temperatures
        y_ref = 22*ones(nz_zones,1); % Zone temperature reference

        %% Prices
        indices_end = min(kk+n_pred-1, length(price_15min_grid)); % Last available price index

        price_grid_h = price_15min_grid(kk:indices_end); % Grid price horizon
        price_dhs_h  = price_15min_dhs(kk:indices_end); % DHS price horizon

        while length(price_grid_h) < n_pred % Pad prices to full prediction horizon
            price_grid_h(end+1) = price_grid_h(end); % Repeat last grid price
            price_dhs_h(end+1)  = price_dhs_h(end); % Repeat last DHS price
        end

        %% Cost
        for k = 1:n_pred % Update linear cost terms for each prediction step

            indices = (k-1)*nx_expand + (1:nx_expand); % Expanded state indices at step k
            f(indices) = -2*C'*Qy*y_ref; % Zone temperature tracking term

            tank_indices = (k-1)*nx_expand + nx_expand; % Tank temperature index at step k
            f(tank_indices) = f(tank_indices) - 2*Q_tank*T_tank_ref; % Tank temperature tracking term

            u_indices = nx_expand*n_pred + (k-1)*nu_expand + (1:nu_expand); % Input indices at step k

            f(u_indices(nu+1)) = f(u_indices(nu+1)) + price_dhs_h(k)*(Ts_hrs/1000); % DHS cost term
            f(u_indices(nu+2)) = f(u_indices(nu+2)) + price_grid_h(k)*(Ts_hrs/(1000*COP)); % Grid cost term

        end

        %% Solve
        model{i}.obj = f; % Update linear objective

        params.OutputFlag = 0; % Suppress Gurobi output
        params.DualReductions = 0; % Disable dual reductions

        result = gurobi(model{i},params); % Solve MPC problem

        if ~isfield(result,'x') % Stop if no solution is returned
            disp(result.status)
            error('Gurobi failed')
        end
        
        %% Extract control
        x_opt = result.x; % Optimal decision vector

        z_end = nx_expand*n_pred; % End index for predicted states
        u_opt = reshape(x_opt(z_end+1:z_end+nu_expand*n_pred),nu_expand,n_pred); % Optimal input sequence

        u_array{i}(:,kk) = u_opt(:,1); % Apply first control move

        u_dhs(i,kk)  = u_array{i}(nu+1,kk); % DHS input [W]
        u_grid(i,kk) = u_array{i}(nu+2,kk); % Grid input [W]
        u_zones{i}(:,kk) = -b*data_{i}.Sq_m_zones'*u_array{i}(1:data_{i}.nu,kk); % Tank temperature decrease from radiator demand
    end

    
    %% STATE UPDATE
    for i = 1:n_buildings % Update state for each building

        data = data_{i}; % Get data for building i
        v_data = data.v_array; % True disturbance data

        u_tmp = u_array{i}(:,kk); % Applied control input at current step

        x_array{i}(:,kk+1) = A_expand{i}*x_array{i}(:,kk) + Bu_expand{i}*u_tmp + Bv_expand{i}*v_data(:,kk); % Simulate next state
    end

end

%% Total supplied power (DHS + Grid) - all buildings

N = n_sim - n_pred;
t = (0:N-1)*Ts_hrs;

Q_dhs_tot  = zeros(1,N);
Q_grid_tot = zeros(1,N);

for i = 1:n_buildings
    nu = data_{i}.nu;

    Q_dhs_tot  = Q_dhs_tot  + u_array{i}(nu+1,1:N); % [W]
    Q_grid_tot = Q_grid_tot + u_array{i}(nu+2,1:N); % [W]
end

if n_days == 1
    t_plot = t;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

figure; hold on;
plot(t_plot, Q_dhs_tot, 'LineWidth',1.5);
plot(t_plot, Q_grid_tot, 'LineWidth',1.5);
plot(t_plot, Q_dhs_tot + Q_grid_tot, 'k','LineWidth',1.5);

legend('DHS','Grid','Total');
title('Total supplied power');
xlabel(xlabel_text);
ylabel('Power [W]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf, fullfile(folder, 'TotalSuppliedPower1.png'));


%% Zone temperatures

N = n_sim - n_pred;

% riktig tidsvektor for states (N+1)
t_x = (0:N) * Ts_hrs;

if n_days == 1
    t_x_plot = t_x;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_x_plot = t_x/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

figure; hold on;

for i = 1:n_buildings
    nz = size(data_{i}.Sq_m_zones,1);
    nx_expand = size(A_expand{i},1);

    C = [eye(nz), zeros(nz, nx_expand-nz)];
    y = C * x_array{i}(:,1:N+1);

    plot(t_x_plot, y', 'LineWidth',1.5)
end

yline(22,'r--') % reference temp
title('Zone temperatures')
xlabel(xlabel_text);
ylabel('Temp  [°C]');
legend('Zone1','Zone2','Zone3','Zone4','Zone5');
grid on;
xlim(x_limits)
ylim([20 30])

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
end

exportgraphics(gcf, fullfile(folder, 'ZoneTemps1.png'));


%% Tank Temps

if n_days == 1
    t_x_plot = t_x;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_x_plot = t_x/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

figure; hold on;
for i = 1:n_buildings
    plot(t_x_plot, x_array{i}(end,1:N+1), 'LineWidth',1.5);
end

legend('T1','T2','T3','T4');
title('Tank temperatures');
xlabel(xlabel_text);
ylabel('Temp  [°C]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
end

exportgraphics(gcf, fullfile(folder, 'TankTemps1.png'));





%% Radiator power

N = n_sim - n_pred;
t_u = (0:N-1)*Ts_hrs;

if n_days == 1
    t_u_plot = t_u;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_u_plot = t_u/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

figure; hold on;
for i = 1:n_buildings
    Q_rad = data_{i}.Sq_m_zones'*u_array{i}(1:data_{i}.nu,1:N);
    plot(t_u_plot, Q_rad, 'LineWidth',1.5);
end

title('Radiator power');
xlabel(xlabel_text);
ylabel('Power [W]');
legend('B1','B2','B3','B4');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf, fullfile(folder, 'RadiatorPowerDHSgrid1.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q_dhs_total  = zeros(1,N);
Q_grid_total = zeros(1,N);

for kk = 1:N
    for i = 1:n_buildings
        nu = data_{i}.nu;

        Q_dhs_total(kk)  = Q_dhs_total(kk)  + u_array{i}(nu+1,kk); % [W]
        Q_grid_total(kk) = Q_grid_total(kk) + u_array{i}(nu+2,kk); % [W]
    end
end

E_dhs_kWh       = sum(Q_dhs_total)* Ts_hrs / 1000;
E_grid_heat_kWh = sum(Q_grid_total) * Ts_hrs / 1000;
E_grid_el_kWh   = E_grid_heat_kWh / COP;

fprintf('Mean DHS power: %.2f kW\n', mean(Q_dhs_total)/1000);
fprintf('Mean grid heat power: %.2f kW\n', mean(Q_grid_total)/1000);
fprintf('Total heat power mean: %.2f kW\n', mean(Q_dhs_total + Q_grid_total)/1000);

fprintf('DHS energy: %.2f kWh\n', E_dhs_kWh);
fprintf('Grid heat energy: %.2f kWh\n', E_grid_heat_kWh);
fprintf('Grid electric energy: %.2f kWh\n', E_grid_el_kWh);

fprintf('Mean DHS price: %.4f NOK/kWh\n', mean(price_15min_dhs(1:N)));
fprintf('Mean grid price: %.4f NOK/kWh\n', mean(price_15min_grid(1:N)));

fprintf('Expected DHS cost: %.2f NOK\n', E_dhs_kWh * mean(price_15min_dhs(1:N)));
fprintf('Expected grid cost approx: %.2f NOK\n', E_grid_el_kWh * mean(price_15min_grid(1:N)));




%% Cumulative cost per building (4 plots)


if n_days == 1
    t_u_plot = t;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_u_plot = t/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

for i = 1:n_buildings

    data = data_{i};
    nu = data.nu;

    % Extract power [W]
    Q_dhs  = u_array{i}(nu+1,1:N);
    Q_grid = u_array{i}(nu+2,1:N);

    % Prices
    price_dhs  = price_15min_dhs(1:N)';
    price_grid = price_15min_grid(1:N)';

    % Cost per step
    cost_dhs  = price_dhs  .* Q_dhs  * Ts_hrs / 1000;
    cost_grid = price_grid .* Q_grid * Ts_hrs / (1000*COP);

    % Cumulative
    cum_dhs   = cumsum(cost_dhs);
    cum_grid  = cumsum(cost_grid);
    cum_total = cum_dhs + cum_grid;

    % Plot
    figure;
    plot(t_u_plot, cum_total, 'LineWidth', 2); hold on;
    plot(t_u_plot, cum_dhs,   'LineWidth', 1.5);
    plot(t_u_plot, cum_grid,  'LineWidth', 1.5);

    grid on;
    xlabel(xlabel_text);
    ylabel('Cost  [NOK]');
    title(['Cumulative cost - Building ', num2str(i)]);
    legend('Total (DHS + Grid)','DHS','Grid');
    xlim(x_limits);

    if n_days == 1
        xticks(0:2:24)
        xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
    end

    % Print final numbers
    fprintf('Building %d:\n', i);
    fprintf('  DHS:   %.2f NOK\n', cum_dhs(end));
    fprintf('  Grid:  %.2f NOK\n', cum_grid(end));
    fprintf('  Total: %.2f NOK\n\n', cum_total(end));

    filename = fullfile(folder, ['CumulativeCosts1_.png', num2str(i), '.png']);
    exportgraphics(gcf, filename);
end


%% Total u_dhs, u_grid and u_zones for all buildings

% allocate storage
u_dhs_total   = zeros(1,N);
u_grid_total  = zeros(1,N);
u_zones_total = zeros(1,N);

for i = 1:n_buildings
    u_dhs_total  = u_dhs_total  + u_dhs(i,1:N);
    u_grid_total = u_grid_total + u_grid(i,1:N);

    % u_zones{i} has one value per room/zone - summing over zones first
    u_zones_total = u_zones_total + sum(u_zones{i}(:,1:N),1);
end

if n_days == 1
    t_u_plot = t;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_u_plot = t/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

figure;
plot(t_u_plot, u_dhs_total,   'LineWidth', 1.5); 
hold on;
plot(t_u_plot, u_grid_total,  'LineWidth', 1.5);
plot(t_u_plot, u_zones_total, 'LineWidth', 1.5);

grid on;
xlabel(xlabel_text);
ylabel('Temp change  [°C/timestep]')
title('Total u\_DHS, u\_grid and u\_zones');
legend('u\_DHS','u\_grid','u\_zones');
xlim(x_limits);

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
end

exportgraphics(gcf, fullfile(folder, 'u_dhs_grid_zones_total1_.png'));

%% u_dhs, u_grid, u_zones per building

if n_days == 1
    t_u_plot = t;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_u_plot = t/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

for i = 1:n_buildings

    % allocating storage
    u_dhs_i   = u_dhs(i,1:N);
    u_grid_i  = u_grid(i,1:N);
    u_zones_i = sum(u_zones{i}(:,1:N),1);

    figure;
    plot(t_u_plot, u_dhs_i,   'LineWidth', 1.5); 
    hold on;
    plot(t_u_plot, u_grid_i,  'LineWidth', 1.5);
    plot(t_u_plot, u_zones_i, 'LineWidth', 1.5);

    grid on;
    xlabel(xlabel_text);
    ylabel('Temp change  [°C/timestep]')
    title(['Building ', num2str(i), ': u\_DHS, u\_grid, u\_zones']);
    legend('u\_DHS','u\_grid','u\_zones');
    xlim(x_limits);

    if n_days == 1
        xticks(0:2:24)
        xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
    end

    filename = fullfile(folder, ['u_dhs_grid_rooms_allbuildings1_', num2str(i), '.png']);
    exportgraphics(gcf, filename);

end

%% Cumulative cost - all buildings

N = n_sim - n_pred;
t = (0:N-1) * Ts_hrs;

Q_dhs_tot  = zeros(1,N);
Q_grid_tot = zeros(1,N);

for i = 1:n_buildings
    nu = data_{i}.nu;

    Q_dhs_tot  = Q_dhs_tot  + u_array{i}(nu+1,1:N); % [W]
    Q_grid_tot = Q_grid_tot + u_array{i}(nu+2,1:N); % [W]
end

if n_days == 1
    t_plot = t;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

price_dhs  = price_15min_dhs(1:N)';
price_grid = price_15min_grid(1:N)';

cost_dhs  = price_dhs  .* Q_dhs_tot  * Ts_hrs / 1000;
cost_grid = price_grid .* Q_grid_tot * Ts_hrs / (1000*COP);

cum_dhs   = cumsum(cost_dhs);
cum_grid  = cumsum(cost_grid);
cum_total = cum_dhs + cum_grid;

figure; hold on;
plot(t_plot, cum_total, 'k', 'LineWidth',2);
plot(t_plot, cum_dhs,   'LineWidth',1.5);
plot(t_plot, cum_grid,  'LineWidth',1.5);

legend('Total','DHS','Grid');
title('Cumulative cost - all buildings');
xlabel(xlabel_text);
ylabel('Cost [NOK]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'})
end

exportgraphics(gcf, fullfile(folder, 'CumulativeCost_AllBuildings1.png'));




%% Price signals

N = n_sim - n_pred;
t_u = (0:N-1) * Ts_hrs;

if n_days == 1
    t_plot = t_u;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t_u/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

price_dhs_plot  = price_15min_dhs(1:N);
price_grid_plot = price_15min_grid(1:N);

figure; hold on;
plot(t_plot, price_dhs_plot,  'LineWidth',1.5);
plot(t_plot, price_grid_plot, 'LineWidth',1.5);

legend('DHS price','Grid price');
title('Energy prices');
xlabel(xlabel_text);
ylabel('NOK/kWh');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf, fullfile(folder, 'Price_Signals1.png'));




%% Total heat extracted from tank per building [kWh]

N = n_sim - n_pred;

fprintf('\nTotal heat extracted from tank per building:\n');

for i = 1:n_buildings

    Q_tank_extract = data_{i}.Sq_m_zones'*u_array{i}(1:data_{i}.nu,1:N);   % [W], summed over zones

    E_tank_extract_kWh = sum(Q_tank_extract) * Ts_hrs / 1000;

    fprintf('Building %d: %.2f kWh\n', i, E_tank_extract_kWh);

end



%% Sum of u_zones per building

fprintf('\nTotal extraction from tank per building:\n');

for i = 1:n_buildings

    total_extract = sum(u_zones{i}(1:N));

    fprintf('Building %d: %.6f °C\n', i, total_extract);

end