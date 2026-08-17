%% LAMBDA IMPLEMENTATION WITH OPTIONAL STRESS TESTS AND PRICE CAP
% Building side MPC


%% Paths
% Add BRCM helper functions.
addpath(genpath('./fcns'));

% Get root path and move to project root.
BRCMRootPath = getBRCMRootPath();
cd(BRCMRootPath);

%% Output folder
folder = 'C:\Users\Malin Storm\Documents\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\Building_All\Lambda\PRESENTATION\PRESENTATION_TEST3';

if ~exist(folder,'dir')
    mkdir(folder)
end

run('DHS_Water_Tank.m') % dependency, init og DHS water tank

%% Building data files
building_data = {'building1_data.mat','building2_data.mat','building3_data.mat','building4_data.mat'};


n_buildings = 4; % Number of buildings
n_days = 1; % Number of simulation days
n_pred = 24;
n_sim = n_days*24*4 + n_pred;  % Total simulation steps plus extra prediction buffer
Ts_hrs = 0.25; % Sampling time [h]


% initial values and params
T_DHS_min = 60;    % [°C]
T_DHS_max = 1000;    % [°C]
Q_DHS_max = 600000; % [W]  (max effect into DHS)
Q_grid    = 300000; % [W]  (max effect from grid)
Qmax = 300000;     % [W]  (max effect bounds) Max from DHS to buildings

% load prices for el grid
price_tbl = readtable('strompriser_30d_NO1.xlsx');
price = price_tbl{:,2}; % choose column 2
price = price(:); % make column vec
price_15min_grid = repelem(price,4); % repeat price 4 times (15 min)
P_reference = mean(price_15min_grid); % monthly based el grid price (nov-dec 30 days)
COP = 3.53; % Coefficient of performance


%% epsilon + lambda init + kappa 

% gradient step size - speed of convergence
kappa = 0.1; %5e-5;%1e-4; %1e-3; %1e-6; % step size for Lambda update
epsilon = 10; %100; %1e-5;  %500%1e2;%1e3; % Convergence tolerance for Q_req - Q_sup
% Initialize lambda using the first grid price
lambda = price_15min_grid(1)/COP*ones(n_buildings,1); 


%% storage allocation per building
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
Q_req_history = zeros(n_buildings,n_sim);

%% BUILD MPC PER BUILDING
for i = 1:n_buildings

    data = load(building_data{i});
    data_{i} = data;

    A = data.A;
    Bu = data.Bu;
    Bv = data.Bv;

    nx = data.nx;
    nu = data.nu;
    nv = data.nv;

    Sq_m_zones = data.Sq_m_zones;
    V = data.tank_size.V;

    Ts = data.Ts_hrs*3600;
    b = Ts/(1000*4186*V);

    A_expand{i} = blkdiag(A,1);
    Bu_expand{i} = [Bu, zeros(nx,2);-b*Sq_m_zones', b, b];
    Bv_expand{i} = [Bv; zeros(1,nv)];

    nx_expand = size(A_expand{i},1);
    nu_expand = size(Bu_expand{i},2);

    n_pred = data.mpc_tuning.n_pred;
    nz_zones = size(Sq_m_zones,1);

    ns_ = nz_zones*n_pred;
    nx_expand_ = nx_expand*n_pred;
    nu_expand_ = nu_expand*n_pred;
    
    %% Prediction matrices
    Omega{i} = getOmega(A_expand{i},n_pred,'lti');
    Psi{i}   = getPsi(A_expand{i},Bu_expand{i},n_pred,'lti');
    Gamma{i} = kron(eye(n_pred),Bv_expand{i});

    %% Output matrices
    I = eye(nz_zones);
    I_nz = zeros(nz_zones,nx_expand);
    I_nz(:,1:nz_zones) = I;
    C = I_nz;
    I_nz_ = kron(eye(n_pred), I_nz);

    %% Bounds 
    T_min = 15; T_max = 30;
    T_min_wt = 35; T_max_wt = 55;
    
    lb_x = -inf*ones(nx_expand_,1);
    ub_x = inf*ones(nx_expand_,1);

    for k = 1:n_pred
        tank_indices = (k-1)*nx_expand + nx_expand;
        lb_x(tank_indices) = T_min_wt;
        ub_x(tank_indices) = T_max_wt;
    end

    u_min = -data.g(1:2:end);
    u_max = data.g(2:2:end);
   
    u_min_expand = [u_min;0;0];
    u_max_expand = [u_max;Qmax;Qmax];

    lb_u = repmat(u_min_expand,n_pred,1);
    ub_u = repmat(u_max_expand,n_pred,1);

    lb_sigma = zeros(ns_,1);
    ub_sigma = inf(ns_,1);

    lb_sigma_tank = zeros(n_pred,1);
    ub_sigma_tank = inf(n_pred,1);

    lb = [lb_x; lb_u; lb_sigma; lb_sigma_tank]; % lower bounds
    ub = [ub_x; ub_u; ub_sigma; ub_sigma_tank]; % upper bounds

    
    %% TEMPERATUR CONSTRAINTS 
    Aineq_temp = [I_nz_, zeros(ns_,nu_expand_), -eye(ns_), zeros(ns_,n_pred);
             -I_nz_, zeros(ns_,nu_expand_), -eye(ns_), zeros(ns_,n_pred)];

    bineq_temp = [T_max*ones(ns_,1);-T_min*ones(ns_,1)];

    %%  ineq CONSTRAINT (TANK + RADIATOR) - with slack

    % allocating storage
    A_tank = zeros(n_pred,length(lb));
    b_tank = zeros(n_pred,1);

    for k = 1:n_pred 
        % placement of dec vars in vector z = [x; u; sigma; sigma_tank]
        u_indices = nx_expand_ + (k-1)*nu_expand + (1:nu);
        tank_indices = (k-1)*nx_expand + nx_expand;
        sigma_tank_indices = nx_expand_ + nu_expand_ + ns_ + k;
        
        % radiator --> tank coupling
        % Sq_m_zones' * u_rad - (1/b)*T_tank - sigma_tank <= -T_min_wt/b
        A_tank(k,u_indices) = Sq_m_zones'; % storing m^2 per zones
        A_tank(k,tank_indices) = -1/b; % storing coefficent
        A_tank(k,sigma_tank_indices) = -1; % storing coefficient
        b_tank(k) = -T_min_wt/b; % right hand side of ineq
    end

    %% INEQUALITIES 
   
    Aineq = [Aineq_temp; A_tank];
    bineq = [bineq_temp; b_tank];

    Aineq_{i} = Aineq;
    bineq_{i} = bineq;

   
    %% Cost weights

    Qy = 1e5; 
    quad_pen = 1e6; %data.mpc_tuning.quad_pen;
    
    T_tank_ref = 40;
    Q_tank = 10;

    E_tank = diag([zeros(nx_expand-1,1);1]);
    
    R = 1e-6*eye(nu_expand);
    R(nu+1,nu+1) = 1e-8;
    R(nu+2,nu+2) = 1e-8;
    
    % Hessian
    H = blkdiag(kron(eye(n_pred), C'*Qy*C + Q_tank*E_tank),kron(eye(n_pred), R),quad_pen*eye(ns_),1e6*eye(n_pred));
    
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

%% Optional price cap storage

% lambda_raw_hist = nan(n_buildings,n_days*96);
% lambda_cap_hist = nan(n_buildings,n_days*96);
% lambda_payment_hist = nan(n_buildings,n_days*96);
% subsidy_price_hist = nan(n_buildings,n_days*96);


%% Optional grid outage stress test

% steps_per_day = 96;
% outage_start = 1*steps_per_day + 1; % Start day 2
% outage_end   = 2*steps_per_day;     % End day 2

% Store original upper bounds for planned outage case
% ub_base = cell(1,n_buildings);
% for i = 1:n_buildings
%     ub_base{i} = model{i}.ub;
% end

% Storage for unexpected outage case
% u_grid_delivered = zeros(n_buildings,n_sim);
% Q_grid_delivered_hist = zeros(1,n_sim);

%% SIMULATION

Dual_Var = 5000; % Maximum number of dual iterations per time step %1500 used when generating results for report;
lambda_dual_hist = nan(n_buildings,Dual_Var,n_days*96); % Store lambda history
kk_plot = 1:96; % Time steps used for plotting
Q_sup_out = zeros(n_buildings,n_sim); % DHS supply delivered to each building

decay = 0.5; % Reduce kappa if oscillations are detected
min_step_size = 1e-9; % Minimum allowed kappa

%% Residual storage
residual_dual_hist = nan(Dual_Var,n_sim); % Dual residual history
residual_building_hist = nan(n_buildings,Dual_Var,n_sim); % Residual per building
iteration_hist = zeros(1,n_sim); % Iterations used at each time step

for kk = 1:(n_sim - n_pred) % Main simulation loop

    % Generate fixed forecast noise for this time step
    w_temp_store = cell(1,n_buildings);
    w_sol_store  = cell(1,n_buildings);

    for ii = 1:n_buildings
        data_tmp = data_{ii};
        v_array_tmp = data_tmp.v_array;

        v_end_tmp = min(kk+n_pred-1,size(v_array_tmp,2));
        v_pred_tmp = v_array_tmp(:,kk:v_end_tmp);

        while size(v_pred_tmp,2) < n_pred
            v_pred_tmp(:,end+1) = v_pred_tmp(:,end);
        end

        w_temp_store{ii} = 0.2*randn(1,size(v_pred_tmp,2));
        w_sol_store{ii}  = 15*randn(1,size(v_pred_tmp,2));
    end

    diff_previous = zeros(n_buildings,1); % Store previous Q_req - Q_sup for oscillation check

    for iteration = 1:Dual_Var % Dual iterations for lambda convergence

        Q_req = zeros(n_buildings,1); % DHS request from each building

        for i = 1:n_buildings % Solve building MPC for current lambda

            data = data_{i}; % Data for building i

            n_pred = data.mpc_tuning.n_pred; % Prediction horizon
            nx_expand = size(A_expand{i},1); % Number of states including tank
            nu_expand = size(Bu_expand{i},2); % Number of inputs including DHS and grid

            nx_expand_ = nx_expand*n_pred; % Total expanded states over horizon
            nu_expand_ = nu_expand*n_pred; % Total expanded inputs over horizon

            %% Predictions
            v_array = data_{i}.v_array; % Disturbance data
            weather_idx = data_{i}.weather_idx; % Weather disturbance indices
            v_end = min(kk + n_pred - 1, size(v_array,2)); % Last disturbance index
            v_pred = v_array(:,kk:v_end); % Disturbance prediction

            while size(v_pred,2) < n_pred % Pad forecast to full horizon
                v_pred(:,end+1) = v_pred(:,end); % Repeat last disturbance value
            end

            w_temp = w_temp_store{i}; % Fixed temperature noise for this time step
            w_sol  = w_sol_store{i}; % Fixed solar noise for this time step

            v_pred(weather_idx.Tamb,:) = v_pred(weather_idx.Tamb,:) + w_temp; % Add temperature noise
            v_pred(weather_idx.SolE,:) = max(0, v_pred(weather_idx.SolE,:) + w_sol); % Add east solar noise
            v_pred(weather_idx.SolS,:) = max(0, v_pred(weather_idx.SolS,:) + w_sol); % Add south solar noise

            v_pred = v_pred(:); % Stack disturbance prediction into one vector

            B_eq = Omega{i}*x_array{i}(:,kk) + Gamma{i}*v_pred; % Dynamic RHS

            nz_zones = size(data.Sq_m_zones,1); % Number of zones
            ns_ = nz_zones*n_pred; % Number of zone slack variables

            n_dec = nx_expand_ + nu_expand_ + ns_ + n_pred; % Number of decision variables
            model{i}.rhs = [bineq_{i}; B_eq]; % Update RHS
            f = zeros(n_dec,1); % Linear cost vector

            I = eye(nz_zones); % Identity matrix for zone temperatures
            I_nz = zeros(nz_zones,nx_expand); % Output selection matrix
            I_nz(:,1:nz_zones) = I; % Select zone temperature states
            C = I_nz; % Output matrix

            y_ref = 22*ones(nz_zones,1); % Zone temperature reference
            
            for k = 1:n_pred % Build cost terms over prediction horizon

                indices = (k-1)*nx_expand + (1:nx_expand); % State indices at step k
                f(indices) = -2*C'*Qy*y_ref; % Zone temperature tracking term

                tank_indices = (k-1)*nx_expand + nx_expand; % Tank temperature index at step k
                f(tank_indices) = f(tank_indices) - 2*Q_tank*T_tank_ref; % Tank temperature tracking term

                u_indices = nx_expand*n_pred + (k-1)*nu_expand + (1:nu_expand); % Input indices at step k

                % Lambda coupling

                price_grid_h = price_15min_grid(kk:kk+n_pred-1); % Grid price prediction horizon

                f(u_indices(nu+1)) = f(u_indices(nu+1)) + lambda(i)*(Ts_hrs/1000); % DHS cost weighted by lambda
                f(u_indices(nu+2)) = f(u_indices(nu+2)) + price_grid_h(k)*(Ts_hrs/(1000*COP)); % Grid electricity cost

            end
            


            %% Optional stress test 1: outage included in the MPC constraints
            
            % The grid input upper bound is set to zero inside the prediction horizon.

            % model{i}.ub = ub_base{i}; % Reset to original input bounds before applying outage
            %
            % for k_pred = 1:n_pred % Loop through prediction horizon
            %     kk_actual = kk + k_pred - 1; % Actual simulation timestep for prediction step k_pred
            %
            %     if kk_actual >= outage_start && kk_actual <= outage_end % Check if step is inside outage period
            %         grid_idx = nx_expand*n_pred + (k_pred-1)*nu_expand + (nu+2); % Grid input index in decision vector
            %         model{i}.ub(grid_idx) = 0; % Disable grid input during outage
            %     end
            % end

            model{i}.rhs = [bineq_{i}; B_eq]; % Update RHS
            model{i}.obj = f; % Update linear objective

            params.OutputFlag = 0; % Suppress Gurobi output
            params.DualReductions = 0; % Disable dual reductions
            result = gurobi(model{i},params); % Solve MPC problem

            x_opt = result.x; % Optimal decision vector

            z_end = nx_expand*n_pred; % End index for predicted states
            u_opt = reshape(x_opt(z_end+1:z_end+nu_expand*n_pred),nu_expand,n_pred); % Optimal input sequence

            u_array{i}(:,kk) = u_opt(:,1); % Apply first control move
            u_dhs(i,kk)  = u_array{i}(nu+1,kk); % DHS request [W]
            u_grid(i,kk) = u_array{i}(nu+2,kk); % Grid input [W]

            format long g % Display full precision

            Q_rad = data_{i}.Sq_m_zones' * u_array{i}(1:data_{i}.nu,kk); % Radiator heat supplied to zones [W]

            Q_req(i) = u_array{i}(nu+1,kk); % Building DHS request [W]
            Q_req_history(i,kk) = u_array{i}(data_{i}.nu+1,kk); % Store DHS request history
           
        end
       
        % Solve DHS MPC and compute supplied heat
        
        %[Q_sup, Q_in] = DHS_MPC_surplus(TankTemp_DHS(kk), Q_req, b_DHS, Qmax, Q_DHS_max,T_DHS_min, T_DHS_max, lambda, kappa);
        [Q_sup, Q_in] = DHS_MPC_final(TankTemp_DHS(kk), Q_req, b_DHS, Qmax, Q_DHS_max,T_DHS_min, T_DHS_max, lambda, kappa);
        
        
      
        lambda_hist(:,kk) = lambda; % Store lambda history
        Q_sup_hist(:,kk)  = Q_sup; % Store DHS supply history
        Q_grid_hist(kk)   = sum(u_grid(:,kk)); % Store total grid usage
        grid_price_hist(kk) = price_15min_grid(kk)/COP; % Effective grid heat price

         
%         % enforce DHS capacity/normalize
%         if sum(Q_sup) > Q_DHS_max
%             Q_sup = Q_sup * (Q_DHS_max / sum(Q_sup));
%         end
% % 
%         LAM_NEW = lambda + kappa*(Q_req-Q_sup); %   eq393
%         
       
        diff = Q_req - Q_sup; % Supply-demand mismatch
   
        residual_dual_hist(iteration, kk) = norm(diff, 2); % Dual residual
        residual_building_hist(:, iteration, kk) = diff; % Residual per building

        if iteration > 1
            sign_change = sign(diff) ~= sign(diff_previous); % Check residual sign changes
            close_to_zero = abs(diff) < 1000; % Ignore large mismatches

            oscillation = any(sign_change & close_to_zero); % Detect oscillation

            if oscillation
                kappa = max(decay*kappa,min_step_size); % Reduce step size if oscillating
            end
        end
        
        diff_previous = diff; % Store residual for next iteration
        % Normal lambda update
        lambda = lambda + kappa*(diff); % Dual update
        lambda = max(lambda,0); % Enforce non-negative lambda

        %% Optional price cap alternatives
        % To activate price cap, comment out the two normal lambda lines above,
        % then uncomment ONE of the alternatives below.

        % Uncapped lambda update
        % lambda_raw = lambda + kappa*diff;
        % lambda_raw = max(lambda_raw,0);

        % Price cap level
        % alpha_consumer_cost = 0.20/COP;
        % lambda_cap = price_15min_grid(kk)/COP + alpha_consumer_cost;
        % lambda_cap = lambda_cap*ones(n_buildings,1);

        % OPTION A: Hard price cap
        % lambda = min(lambda_raw,lambda_cap);
        % lambda_payment = lambda;
        % subsidy_price = max(0,lambda_raw - lambda_cap);

        % OPTION B: Payment cap only
        % lambda_payment = min(lambda_raw,lambda_cap);
        % subsidy_price = max(0,lambda_raw - lambda_cap);
        % lambda = lambda_raw;

        % Store price cap values
        % lambda_raw_hist(:,kk) = lambda_raw;
        % lambda_cap_hist(:,kk) = lambda_cap;
        % lambda_payment_hist(:,kk) = lambda_payment;
        % subsidy_price_hist(:,kk) = subsidy_price;
        
        lambda_dual_hist(:,iteration,kk) = lambda; % Store lambda iteration history
        
        Q_sup_out(:,kk) = Q_sup; % DHS supply delivered to buildings


        if norm(Q_sup - Q_req,2) < epsilon % Stop dual iterations if converged
            break;
        end
                
    end
    
    iteration_hist(kk) = iteration; % Store number of dual iterations 

    for i = 1:n_buildings

        v_array = data_{i}.v_array; % True disturbance data

        u_tmp = u_array{i}(:,kk); % Building MPC input
        u_tmp = u_tmp(:); % Convert to column vector

        u_tmp(nu+1) = Q_sup_out(i,kk); % Replace DHS request with supplied DHS heat
        
        %% Optional stress test 2: unexpected grid outage

        % if kk >= outage_start && kk <= outage_end
        %     u_grid_delivered(i,kk) = 0;
        %     u_tmp(nu+2) = 0; % Grid unavailable in real system
        % else
        %     u_grid_delivered(i,kk) = u_grid(i,kk);
        %     u_tmp(nu+2) = u_grid(i,kk);
        % end

        x_array{i}(:,kk+1) = A_expand{i}*x_array{i}(:,kk) + Bu_expand{i}*u_tmp + Bv_expand{i}*v_array(:,kk); % Update building state

    end
    
    % TANK UPDATE EQUATION
    %TankTemp_DHS(kk+1) = DHSTank(TankTemp_DHS(kk), Q_sup_out(:,kk), b_DHS, T_DHS_min, T_DHS_max);
    TankTemp_DHS(kk+1) = DHSTank(TankTemp_DHS(kk), Q_sup_out(:,kk), Q_in,b_DHS, T_DHS_min, T_DHS_max);
    % Q_grid_delivered_hist(kk) = sum(u_grid_delivered(:,kk));
end
    
% save simulation data
save(fullfile(folder,'FULL_workspace_dist_RESIDUAL.mat'));


kk_show = 20;

figure; hold on;

for i = 1:n_buildings
    r_i = squeeze(residual_building_hist(i,:,kk_show));
    valid_iter = ~isnan(r_i);
    plot(find(valid_iter), r_i(valid_iter), 'LineWidth', 1.5);
end

yline(0, 'k--');

legend('Building 1','Building 2','Building 3','Building 4');
title(['Q_{req} - Q_{sup} per dual iteration, kk = ', num2str(kk_show)]);
xlabel('Dual iteration');
ylabel('Q_{req} - Q_{sup} [W]');
grid on;

exportgraphics(gcf, fullfile(folder, ['DualResidual_PerBuilding_kk', num2str(kk_show), '.png']));


%% Plot lambda convergence for kk = m:n

kk_plot = 19:20;

figure; hold on;

for i = 1:n_buildings

    lambda_tmp = [];

    for j = 1:length(kk_plot)

        kk_idx = kk_plot(j);

        valid_idx = lambda_dual_hist(i,:,kk_idx) ~= 0;

        lambda_piece = squeeze(lambda_dual_hist(i,valid_idx,kk_idx));
        lambda_piece = lambda_piece(:)';

        lambda_tmp = [lambda_tmp,lambda_piece];

    end

    plot(lambda_tmp,'LineWidth',1.5);

end

legend('Building 1','Building 2','Building 3','Building 4');

title('Convergence of \Lambda, kk = 19:20');

xlabel('Total dual iteration');

ylabel('\lambda [NOK/kWh]');

grid on;

exportgraphics(gcf, fullfile(folder, 'Lambda_Convergence_3.png'));

% %% Lambda vs grid price
N = n_sim - n_pred;
t = (0:N-1) * 0.25;

if n_days == 1
    t_plot = t;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

N = min([n_sim - n_pred, size(lambda_hist,2), length(grid_price_hist)]);

lambda_mean = mean(lambda_hist(:,1:N),1);
grid_price = price_15min_grid(1:N)/COP;

% Old fixed DHS price
lambda_fixed = 1.1723*ones(1,N);

figure; hold on;

plot(t_plot(1:N), lambda_mean, 'LineWidth',1.5);
plot(t_plot(1:N), grid_price, 'LineWidth',1.5);
plot(t_plot(1:N), lambda_fixed, '--', 'LineWidth',1.5);

legend('Mean \lambda','Grid price / COP','Fixed DHS price = 1.2');
title('\lambda compared with grid price');
xlabel(xlabel_text);
ylabel('Price [NOK/kWh]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24)
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf, fullfile(folder, 'Lambda_vs_Grid7.png'));
% 
% %%%%%%%%%%%%%%%% Total Request %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

% plot(t_plot, sum(u_dhs(:,1:N),1),'LineWidth',1.5);
% plot(t_plot, sum(u_grid(:,1:N),1),'LineWidth',1.5);
% plot(t_plot, sum(u_dhs(:,1:N) + u_grid(:,1:N),1),'k','LineWidth',1.5);
Q_dhs_plot  = sum(u_dhs(:,1:N),1);
Q_grid_plot = sum(u_grid(:,1:N),1);

plot(t_plot(1:N), Q_dhs_plot(1:N), 'LineWidth',1.5);
plot(t_plot(1:N), Q_grid_plot(1:N), 'LineWidth',1.5);
plot(t_plot(1:N), Q_dhs_plot(1:N) + Q_grid_plot(1:N),'k','LineWidth',1.5);

legend('DHS request','Grid','Total request');
title('Total requested power');
xlabel(xlabel_text);
ylabel('Power [W]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end


exportgraphics(gcf, fullfile(folder, 'TotalRequest7.png'));
% 
% 
% %%%%%%%%%%%%% Total supply %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_dhs = sum(Q_sup_out,1);
Q_grid = sum(u_grid,1);

N2 = min(length(Q_dhs), length(Q_grid));

figure; hold on;

plot(t_plot(1:N2), Q_dhs(1:N2),'LineWidth',1.5);
plot(t_plot(1:N2), Q_grid(1:N2),'LineWidth',1.5);
plot(t_plot(1:N2), Q_dhs(1:N2) + Q_grid(1:N2),'k','LineWidth',1.5);

legend('DHS supply','Grid','Total supply');
title('Total supplied power');
xlabel(xlabel_text);
ylabel('Power [W]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end
exportgraphics(gcf, fullfile(folder, 'TotalSupply7.png'));

% 
% %%%%%%%%%%%%%%%%%%% Qreq - Qsup TOTAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_req_tot = sum(u_dhs,1);
Q_sup_tot = sum(Q_sup_out,1);

N = min(length(Q_req_tot), length(Q_sup_tot));

figure;
plot(t_plot(1:N), Q_req_tot(1:N) - Q_sup_tot(1:N),'LineWidth',1.5);

yline(0,'k--');

title('Q_{req} - Q_{sup}');
xlabel(xlabel_text);
ylabel('Power [W]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf, fullfile(folder, 'Qreq_minus_QSup7.png'));
% 
% %%%%%%%%% Qreq-Qsup per building %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure; hold on;

for i = 1:n_buildings
    N_i = min(length(u_dhs(i,:)), length(Q_sup_out(i,:)));

    plot(t_plot(1:N_i),u_dhs(i,1:N_i) - Q_sup_out(i,1:N_i),'LineWidth',1.5);
end

yline(0,'k--');

title('Q_{req} - Q_{sup} per building');
xlabel(xlabel_text);
ylabel('Power [W]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf, fullfile(folder, 'Qreq_minus_Qsup_perBuilding7.png'));
% 
% 
% 
% %%%%%%%%%%%%%%%%%%% Zone Temps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

for i = 1:n_buildings
    nz = size(data_{i}.Sq_m_zones,1);
    nx_expand = size(A_expand{i},1);

    C = [eye(nz),zeros(nz,nx_expand-nz)];
    y = C*x_array{i};

    N_i = min([length(t_plot),size(y,2)]);

    plot(t_plot(1:N_i),y(:,1:N_i)','LineWidth',1.5);
end

yline(22,'r--');

title('Zone temperatures');
xlabel(xlabel_text);
ylabel('Temp [°C]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf,fullfile(folder,'Zone_Temperatures7.png'));
% 
% %%%%%%%%%%%%%%%%%%%%%% TANK TEMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

N_tank = min([length(t_plot),length(TankTemp_DHS)]);

for i = 1:n_buildings
    N_i = min([N_tank,size(x_array{i},2)]);
    plot(t_plot(1:N_i),x_array{i}(end,1:N_i),'LineWidth',1.5);
end

plot(t_plot(1:N_tank),TankTemp_DHS(1:N_tank),'k','LineWidth',1.5);

legend('T1','T2','T3','T4','DHS');
title('Tank temperatures');
xlabel(xlabel_text);
ylabel('Temp [°C]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

disp(isgraphics(gcf));
drawnow;
exportgraphics(gcf,fullfile(folder,'TankTemp7.png'));
% 
% %%%%%%%%%%%%%%%%%% RADIATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = n_sim-n_pred;
t_u = (0:N-1)*Ts_hrs;

if n_days == 1
    t_plot = t_u;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t_u/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

figure; hold on;

for i = 1:n_buildings
    Q_rad = data_{i}.Sq_m_zones'*u_array{i}(1:data_{i}.nu,:);

    N_i = min([length(t_plot),length(Q_rad)]);

    plot(t_plot(1:N_i),Q_rad(1:N_i),'LineWidth',1.5);
end

title('Radiator power');
xlabel(xlabel_text);
ylabel('Power [W]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf,fullfile(folder,'TotalRadiatorPower7.png'));

% %%%%%%%%%%%5 Cum Cost per building %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
t_u = (0:N-1)*Ts_hrs;

if n_days == 1
    t_plot = t_u;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t_u/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

for i = 1:n_buildings

    %cost_dhs_i  = lambda_hist(i,1:N) .* Q_sup_hist(i,1:N)*Ts_hrs/1000;
    N_i = min([N,size(lambda_hist,2),size(Q_sup_hist,2),size(u_grid,2),length(grid_price_hist),length(t_plot)]);

    cost_dhs_i = lambda_hist(i,1:N_i).*Q_sup_hist(i,1:N_i)*Ts_hrs/1000;

    %cost_grid_i = (grid_price_hist(1:N)./COP) .* Q_grid_hist(i,1:N)*Ts_hrs / 1000;
    cost_grid_i = grid_price_hist(1:N_i).*u_grid(i,1:N_i)*Ts_hrs/1000;

    cum_dhs_i = cumsum(cost_dhs_i);
    cum_grid_i = cumsum(cost_grid_i);
    cum_tot_i = cumsum(cost_dhs_i + cost_grid_i);

    figure; hold on;
%     plot(t_plot, cum_dhs_i,'LineWidth',1.5);
%     plot(t_plot, cum_grid_i,'LineWidth',1.5);
%     plot(t_plot, cum_tot_i,'k','LineWidth',1.5);

    plot(t_plot(1:N_i),cum_dhs_i,'LineWidth',1.5);
    plot(t_plot(1:N_i),cum_grid_i,'LineWidth',1.5);
    plot(t_plot(1:N_i),cum_tot_i,'k','LineWidth',1.5);

    legend('DHS','Grid','Total');
    title(['Cumulative cost - Building ',num2str(i)]);
    xlabel(xlabel_text);
    ylabel('Cost [NOK]');
    grid on;
    xlim(x_limits);

    if n_days == 1
        xticks(0:2:24);
        xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
    end

    exportgraphics(gcf,fullfile(folder,['Cumulative_Cost_Building7_',num2str(i),'.png']));
end

% %% Cumulative Cost plotting

N = n_sim-n_pred;
t_u = (0:N-1)*Ts_hrs;

if n_days == 1
    t_plot = t_u;
    xlabel_text = 'Time [hours]';
    x_limits = [0 24];
else
    t_plot = t_u/24;
    xlabel_text = 'Time [days]';
    x_limits = [0 n_days];
end

% Q_dhs_tot  = sum(Q_sup_hist(:,1:N),1);
% Q_grid_tot = sum(Q_grid_hist(:,1:N),1);

N_tot = min([N,size(Q_sup_hist,2),size(u_grid,2),size(lambda_hist,2),length(grid_price_hist),length(t_plot)]);

Q_dhs_tot = sum(Q_sup_hist(:,1:N_tot),1);
Q_grid_tot = sum(u_grid(:,1:N_tot),1);
%cost_dhs  = sum(lambda_hist(:,1:N) .* Q_sup_hist(:,1:N),1)*Ts_hrs/1000;
cost_dhs = sum(lambda_hist(:,1:N_tot).*Q_sup_hist(:,1:N_tot),1)*Ts_hrs/1000;
%cost_grid = (grid_price_hist(1:N)./COP) .* Q_grid_tot*Ts_hrs/1000;
cost_grid = grid_price_hist(1:N_tot).*Q_grid_tot*Ts_hrs/1000;

cum_dhs = cumsum(cost_dhs);
cum_grid = cumsum(cost_grid);
cum_tot = cumsum(cost_dhs + cost_grid);

figure; hold on;
% plot(t_plot, cum_dhs,'LineWidth',1.5);
% plot(t_plot, cum_grid,'LineWidth',1.5);
% plot(t_plot, cum_tot,'k','LineWidth',1.5);

plot(t_plot(1:N_tot),cum_dhs,'LineWidth',1.5);
plot(t_plot(1:N_tot),cum_grid,'LineWidth',1.5);
plot(t_plot(1:N_tot),cum_tot,'k','LineWidth',1.5);

legend('DHS','Grid','Total');
title('Cumulative cost');
xlabel(xlabel_text);
ylabel('Cost [NOK]');
grid on;
xlim(x_limits);

if n_days == 1
    xticks(0:2:24);
    xticklabels({'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});
end

exportgraphics(gcf,fullfile(folder,'Cumulative_Cost7.png'));


%% Optional price cap plot: raw lambda, payment lambda, cap and grid

% N = n_days*96;
% t = (0:N-1)*Ts_hrs;
%
% if n_days == 1
%     t_plot = t;
%     xlabel_text = 'Time [hours]';
%     x_limits = [0 24];
% else
%     t_plot = t/24;
%     xlabel_text = 'Time [days]';
%     x_limits = [0 n_days];
% end
%
% lambda_raw_mean = mean(lambda_raw_hist(:,1:N),1,'omitnan');
% lambda_payment_mean = mean(lambda_payment_hist(:,1:N),1,'omitnan');
% lambda_cap_mean = mean(lambda_cap_hist(:,1:N),1,'omitnan');
% grid_price = price_15min_grid(1:N)/COP;
%
% figure; hold on;
% plot(t_plot, lambda_raw_mean, 'LineWidth',1.5);
% plot(t_plot, lambda_payment_mean, 'LineWidth',1.5);
% plot(t_plot, lambda_cap_mean, '--', 'LineWidth',1.5);
% plot(t_plot, grid_price, 'LineWidth',1.5);
% legend('\lambda raw','\lambda payment','Price cap','Grid price / COP');
% title('\lambda raw, payment price, cap and grid price');
% xlabel(xlabel_text);
% ylabel('Price [NOK/kWh]');
% grid on;
% xlim(x_limits);
% exportgraphics(gcf, fullfile(folder, 'Lambda_raw_payment_cap_grid.png'));


%% Optional price cap plot: cumulative payment and subsidy

% cost_DHS_customer_t = sum(lambda_payment_hist(:,1:N).*Q_sup_hist(:,1:N),1)*Ts_hrs/1000;
% cost_DHS_raw_t = sum(lambda_raw_hist(:,1:N).*Q_sup_hist(:,1:N),1)*Ts_hrs/1000;
% subsidy_t = sum(subsidy_price_hist(:,1:N).*Q_sup_hist(:,1:N),1)*Ts_hrs/1000;
%
% cum_DHS_customer = cumsum(cost_DHS_customer_t);
% cum_DHS_raw = cumsum(cost_DHS_raw_t);
% cum_subsidy = cumsum(subsidy_t);
%
% figure; hold on;
% plot(t_plot, cum_DHS_customer, 'LineWidth',1.5);
% plot(t_plot, cum_subsidy, 'LineWidth',1.5);
% plot(t_plot, cum_DHS_raw, 'LineWidth',1.5);
% legend('Customer DHS payment','State subsidy','DHS raw value');
% title('Cumulative payment and subsidy with price cap');
% xlabel(xlabel_text);
% ylabel('Cost [NOK]');
% grid on;
% xlim(x_limits);
% exportgraphics(gcf, fullfile(folder, 'Cumulative_payment_subsidy_pricecap.png'));

%% Optional stress test B plot: requested vs delivered grid

% figure; hold on;
% plot(t_plot, sum(u_grid(:,1:N),1), 'LineWidth',1.5);
% plot(t_plot, Q_grid_delivered_hist(1:N), '--', 'LineWidth',1.5);
% legend('Requested grid','Delivered grid');
% title('Requested and delivered grid power');
% xlabel(xlabel_text);
% ylabel('Power [W]');
% grid on;
% xlim(x_limits);
% exportgraphics(gcf, fullfile(folder, 'Requested_vs_Delivered_Grid.png'));

% Parts of the plotting code, alternative lambda-based solution approaches, comments, outage scenarios and price-cap scenarios were developed
% with assistance from ChatGPT and subsequently reviewed, modified and validated by the author.
