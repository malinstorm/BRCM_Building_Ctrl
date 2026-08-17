%% CALCULATING HEAT DEMAND FOR THE BUILDINGS

clear; clc; close all; % remove variables from WS, clear CMD window, closes all open firgures

%% Paths
% Add all folders inside ./fcns to the MATLAB path.
% fcns contains functions such as getOmega, getPsi,...
addpath(genpath('./fcns'));

% Get root path of the BRCM toolbox/project.
BRCMRootPath = getBRCMRootPath();

% Navigate to folder
cd(BRCMRootPath);
folder = 'C:\Users\Malin Storm\Documents\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\Buildings\3\1_day';

% Global debug level used by BRCM
global g_debugLvl
g_debugLvl = 0;

%% Building
demoBuilding = 'DemoBuilding'; % name of building model

% paths to Thermal data and External Heat Flux data
thermalModelDataDir = [BRCMRootPath,filesep,'BuildingData',filesep,demoBuilding,filesep,'ThermalModel'];
EHFModelDataDir     = [BRCMRootPath,filesep,'BuildingData',filesep,demoBuilding,filesep,'EHFM'];

% Create object and load thermal data for object
B = Building(demoBuilding);
B.loadThermalModelData(thermalModelDataDir);

%% Declare ALL EHF used by Building Model
B.declareEHFModel('BuildingHull.m',[EHFModelDataDir,filesep,'buildinghull'],'BuildingHull');
%B.declareEHFModel('AHU.m',          [EHFModelDataDir,filesep,'ahu'],'AHU1');
B.declareEHFModel('InternalGains.m',[EHFModelDataDir,filesep,'internalgains'],'IG');
%B.declareEHFModel('BEHeatfluxes.m', [EHFModelDataDir,filesep,'BEHeatfluxes'],'TABS');
B.declareEHFModel('Radiators.m',    [EHFModelDataDir,filesep,'radiators'],'Rad');

%% Generate model
% Continuous-time thermal model
B.generateThermalModel;
% Complete building building model including EHF model
B.generateBuildingModel;

%% Discretization
Ts_hrs = 0.25; % Sampling time
B.building_model.setDiscretizationStep(Ts_hrs); % Set discretization step (step size)
B.building_model.discretize(); % discretize

%% System matrices
A = B.building_model.discrete_time_model.A; % Discrete-time state matrix
Bu = B.building_model.discrete_time_model.Bu; % Discrete-time input matrix for controllable inputs
Bv = B.building_model.discrete_time_model.Bv; % Discrete-time disturbance matrix
   

%% Dimension and identifiers
% Print identifiers for states, inputs, and disturbances.
% Checking the names and ordering of variables.
B.building_model.printIdentifiers;


%% Weather (BuildingHull.m)
% Disturbance indices for ambient temperature and solar radiation.
% These indices are used later to fill the disturbance vector v_array.
indices_Tamb_v = getIdIndex('v_Tamb',B.building_model.identifiers.v);
indices_SolE_v = getIdIndex('v_solGlobFac_E',B.building_model.identifiers.v);
indices_SolN_v = getIdIndex('v_solGlobFac_N',B.building_model.identifiers.v);
indices_SolS_v = getIdIndex('v_solGlobFac_S',B.building_model.identifiers.v);

%% Internal Gains 

% Disturbance indices for internal gains in each zone.
% These are used to assign heat gains from people to the correct disturbance channels.
% indices_IG_zone1 = getIdIndex('v_IG_ZoneOne',B.building_model.identifiers.v);
% indices_IG_zone2 = getIdIndex('v_IG_ZoneTwo',B.building_model.identifiers.v);
% indices_IG_zone3 = getIdIndex('v_IG_ZoneThree',B.building_model.identifiers.v);
% indices_IG_zone4 = getIdIndex('v_IG_ZoneFour',B.building_model.identifiers.v);
% indices_IG_zone5 = getIdIndex('v_IG_ZoneFive',B.building_model.identifiers.v);



%% Radiator bounds from BRCM

% HEAT FLUX - min/max W/m^2 control input
% The allowed radiator heat flux in each zone (5).
constraintsParameters.Rad.Q_rad_Radiator1_min = 0;
constraintsParameters.Rad.Q_rad_Radiator1_max = 150; 

constraintsParameters.Rad.Q_rad_Radiator2_min = 0;
constraintsParameters.Rad.Q_rad_Radiator2_max = 150;

constraintsParameters.Rad.Q_rad_Radiator3_min = 0;
constraintsParameters.Rad.Q_rad_Radiator3_max = 150;

constraintsParameters.Rad.Q_rad_Radiator4_min = 0;
constraintsParameters.Rad.Q_rad_Radiator4_max = 150;

constraintsParameters.Rad.Q_rad_Radiator5_min = 0;
constraintsParameters.Rad.Q_rad_Radiator5_max = 150;


%% Linear constraint ineq matrices (Fx = state constRaint + Fu input constraint <= g)
% Get constraint matrices from the BRCM building model.
% Only g is used later to extract input bounds.
[Fx,Fu,~,g] = B.building_model.getConstraintsMatrices(constraintsParameters);



%% Simulation and prediction parameters

% N simulation days
n_days = 1;

% N simulation steps
% 24 hours/day and 4 samples/hour <-- Ts = 0.25 h.
n_sim = n_days*24*4; % n_days*24 hrs*4 times/hour

% N thermal zones in the building.
nz_zones = 5;

% Prediction horizon for MPC
% n_pred = 24 means 24 prediction steps = 6 hours with Ts = 0.25 h.
n_pred = 24; 


%% Model dimensions
nx = size(A,1); % N thermal states N x N states 
nu = size(Bu,2); % Number of ctrl inputs --> 5 ctrl variables
ny = nz_zones; % Number of measured zone temepratures
nv = size(Bv, 2); % N disturbance inputs --> 10 dist variables

%% Prediction horizon
nx_ = nx*n_pred; % Total pred states over horizon
nu_ = nu*n_pred; % Total ctrl variables over horizon
ns_ = nz_zones*n_pred; % Total zone temperatures over horizon

%% Square meter zones
Sq_m_zones = zeros(nu,1); % allocate storage for square meter per zone
for i = 1:nu
    identifier = sprintf('Z%04d',i); % column for zone values
    Sq_m_zones(i) = B.thermal_model_data.getValue(identifier,'area'); 
end

%% WEATHER DATA
% allocate storage - solar radiation & ambient temperatures
SolRad_E = zeros(1,n_sim);
SolRad_S = zeros(1,n_sim);
SolRad_N = zeros(1,n_sim);

T_amb = zeros(1,n_sim);

for i = 1:n_sim

   % Synthetic solar radiation profile
   % The sine function creates a daily pattern
   % max(0, ...) makes sure radiation is never negative
   sine = max(0,(600*sin(2*pi*i/(24*4) + 3*pi/2)-100));

   % emulate diffrent angles of Sunlight due to Earth revolution
   % East-facing solar radiation.
   SolRad_E(i) = sine*0.7;  
   % South-facing solar radiation.
   SolRad_S(i) = sine;
   % North-facing solar radiation is set to zero.
   SolRad_N(i) = 0; % azimuth

   % Synthetic ambient temperature profile.
   % Temperature varies during the day using a sinusoidal function.
   T_amb(i) = -2.5 + 7.5*sin(2*pi*i/(24*4)+3*pi/2); 
end


%% Internal Gains
%Internal_Gains = 5*ones(1,n_sim); % constant ig 5 W/m^2


% load('ihg_model.mat'); % load data from internal gains model/ occypancy profile
% 
% p_array = internalHeatGain.p_array; % occupancy profile from igh.mat
% p_array = p_array(1:nz_zones,:); % Keep only the occupancy profiles for the modelled zones (5 rows, 96 timesteps (n_days = 1))
% 
% p_array = repmat(p_array,1,ceil(n_sim/96)); % repat matrix m x n times 
% p_array = p_array(:,1:n_sim); % ensure only comlumns needed for simulation (for n_days = 1.5, etc)
% 
% Q_IG = 100*p_array; % 100 W scaling factor - scale occupancy profile to Internal Heat Gain
% q_internalgain_m2 = zeros(nu,n_sim); % allocate storage for W/m^2
% 
% for i = 1:nu
%    q_internalgain_m2(i,:) = Q_IG(i,:)./Sq_m_zones(i); % divide by sq m --> [W/m^2]
% end


%% Min/max radiator input bounds [W/m^2]

u_min = -g(1:2:end); % Extract lower radiator bounds from BRCM constraint vector
u_max =  g(2:2:end); % Extract upper radiator bounds from BRCM constraint vector

% Maximum heating power available from (DHS) and grid [W]
%Qmax_DHS = 80000;
%Qmax_grid = 80000;

% Extend input bounds with DHS and grid power variables
%u_min_expand = [u_min; 0; 0]; 
%u_max_expand = [u_max; Qmax_DHS; Qmax_grid];

%% Zone temperature limits [Celsius]
% Min/max allowed zone temperature
T_min = 15;
T_max = 30;


%% Prediction matrices
% The model is linear time-invariant (LTI), meaning A, Bu and Bv are constant.
% This allows the predicted state trajectory over the horizon to be written as:
%
% X = Omega*x0 + Psi*u + Gamma*v
%
% instead of propagating the dynamics step by step inside the optimizer.

Omega = getOmega(A,n_pred,'lti'); % Maps current state x0 to predicted states X
Psi   = getPsi(A,Bu,n_pred,'lti'); % Maps future control inputs u to predicted states X
Gamma = kron(eye(n_pred),Bv); % Maps future disturbances v to predicted states X

%% Cost weights/parameters

% Weight on zone temperature tracking error.
% High value means the MPC strongly tries to keep temperatures near y_ref
Qy = 1000; % high weight put on output error - penalty on deviation from reference

% Small penalty on all thermal states.
Qx = 1e-4*eye(nx);

% Penalty on control input.
% Small R means heat usage is not strongly penalized - heat is cheap
R = 1e-6*eye(nu);

% quadratic penalty - help increase convexity of the quadratic part of the cost function --> increases curvature
% This variable is not directly used in the current Hessian definition.
quad_pen = 50;

%% Initial 

% Initial state vector - All thermal states are initialized to 22 degrees.
x0 = 22*ones(nx,1);


% Zone temperature reference/setpoint.
y_ref = 22*ones(nz_zones,1); % setpoint for the 5 zones

% Lower and upper bounds for predicted states.
% These are loose bounds, so the real comfort limits are handled by Aineq_temp.
lb_x = -1000*ones(nx*n_pred,1);
ub_x = 1000*ones(nx*n_pred,1);



% Lower and upper bounds for future control inputs.
lb_u = repmat(u_min,n_pred,1);
ub_u = repmat(u_max,n_pred,1);


% Combined lower and upper bounds for optimization decision variables.
% Decision variable ordering: [predicted states; predicted controls].
lb = [lb_x;lb_u];
ub = [ub_x;ub_u];

%% Slacks & constraints

% Aineq * x <= bineq
% l(xi) = max{0, Tmin_i - xi, xi - Tmax_i}

% Aineq * x <= bineq -->
% [I 0 -I;-I 0 -I;0 0 -I][x; u; sigma] <= [Tmax; -Tmin; 0]
% [I 0 -I;-I 0 -I[x; u; sigma] <= [Tmax; -Tmin]
% sigma >= 0 --> lb_sigma
% sigma >= 0

%% Output selection matrix
I = eye(nz_zones); % identity num zones matrix (5x5 identity)

% Selects the first 5 states, the zone temperatures.
I_nz = zeros(nz_zones,nx);  
I_nz(:,1:nz_zones) = I; % 5 first columns = I, extracts the 5 zone temperatures from the N-state vector

% Output matrix.
% Zone temperatures are calculated as y = C*x.
C = I_nz; % output measure matrix --> y = Cx

% Block-diagonal version of C over the full prediction horizon.
I_nz_ = kron(eye(n_pred), I_nz); % kron matrix selected temperatures


% Temperature inequality constraints:
% I_nz_ * X <= T_max
% -I_nz_ * X <= -T_min
% Together these enforce:
% T_min <= zone temperatures <= T_max
% Aineq_temp * x <= bineq
Aineq_temp = [I_nz_,zeros(ns_,nu_);-I_nz_,zeros(ns_,nu_)]; % Aineq upper/lower constraints
bineq = [T_max*ones(ns_,1);-T_min*ones(ns_,1)];


% Equality constraint for Gurobi. Prediction model: X = Omega*x0 + Psi*U + Gamma*V
% Rewritten as: X - Psi*U = Omega*x0 + Gamma*V
% Decision variables are [X; U], so A_eq*[X;U] = X - Psi*U.
A_eq = [eye(nx_), -Psi];
% Initial RHS. Updated in MPC loop with current state and predicted disturbances.
B_eq = Omega*x0;


% Hessian matrix for the quadratic cost function.
% The state block penalizes zone temperature tracking and weakly regularizes all states.
% The input block penalizes radiator control effort.
% Consists of a state block Hx = (Qx + C'*Qy*C) and input block Hy = R

H = blkdiag(kron(eye(n_pred),Qx + C'*Qy*C),kron(eye(n_pred),R));


%% Allocate storage

nx_f = size(H,1); % Number of optimization variables

f = zeros(nx_f,1); % Linear cost vector for Gurobi

x_array = zeros(nx,n_sim+1); % State trajectory over simulation
y_array = zeros(ny,n_sim+1); % Zone temperature trajectory
u_array = zeros(nu,n_sim);   % Applied radiator inputs
v_array = zeros(nv,n_sim);   % Disturbance time series


% initial conditions loaded into storage matrices
x_array(:,1) = x0;
y_array(:,1) = C*x0;

% Extract optimal sequences - for z and u - (end of respective part of vector)
z_end = nx*n_pred; 
u_end = z_end + nu*n_pred;

%% Loading values into V_array

% Assign weather disturbances to correct disturbance channels.
v_array(indices_Tamb_v,:) = T_amb;
v_array(indices_SolE_v,:) = SolRad_E;
v_array(indices_SolS_v,:) = SolRad_S;
v_array(indices_SolN_v,:) = SolRad_N;

% Old constant internal gain alternative
%v_array(indices_IG_zone1,:) = Internal_Gains;
%v_array(indices_IG_zone2,:) = Internal_Gains;
%v_array(indices_IG_zone3,:) = Internal_Gains;
%v_array(indices_IG_zone4,:) = Internal_Gains;
%v_array(indices_IG_zone5,:) = Internal_Gains;

% Assign internal gains from occupancy model to each zone
% v_array(indices_IG_zone1,:) = q_internalgain_m2(1,:);
% v_array(indices_IG_zone2,:) = q_internalgain_m2(2,:);
% v_array(indices_IG_zone3,:) = q_internalgain_m2(3,:);
% v_array(indices_IG_zone4,:) = q_internalgain_m2(4,:);
% v_array(indices_IG_zone5,:) = q_internalgain_m2(5,:);

% Initial disturbance vector.
v0 = v_array(:,1);

%% DEFINE GUROBI
% Quadratic cost matrix for Gurobi.
model.Q = sparse(H);

% Constraint matrix for Gurobi.
model.A = sparse([Aineq_temp;A_eq]);% total restriction matrix 

% Constraint senses:
% '<' for temperature inequality constraints, '=' for equality constraints
model.sense =[repmat('<',size(Aineq_temp,1),1);repmat('=',size(A_eq,1),1)];

% Lower and upper bounds for decision variables.
model.lb = lb; % set the lower bound (lb <= x <= ub)
model.ub = ub; % set the upper bound


%% MPC LOOP
for kk = 1:n_sim

    v_end = min(kk + n_pred - 1, n_sim); % Last disturbance index within simulation length
    v_pred = v_array(:,kk:v_end); % Disturbance forecast over available horizon

    while size(v_pred,2) < n_pred % Pad forecast to full prediction horizon
        v_pred(:,end+1) = v_pred(:,end); % Repeat last disturbance value
    end

    % generating random forecast noise each iteration
    w_temp = 0.2*randn(1,size(v_pred,2));   % ambient temperature noise
    w_sol  = 15*randn(1,size(v_pred,2));    % solar radiation noise

    % add forecast noise
    v_pred(indices_Tamb_v,:) = v_pred(indices_Tamb_v,:) + w_temp;
    v_pred(indices_SolE_v,:) = max(0, v_pred(indices_SolE_v,:) + w_sol);
    v_pred(indices_SolS_v,:) = max(0, v_pred(indices_SolS_v,:) + w_sol);

    % vectorize after all edits
    v_pred = v_pred(:); % Stack forecast disturbances into one vector

    % Equality constraint update
    B_eq = Omega*x_array(:,kk) + Gamma*v_pred; % dynamic RHS

    % Reset linear cost term
    f(:) = 0;

    for k = 1:n_pred % Add tracking term for each prediction step
        indices = (k-1)*nx + (1:nx); % State indices at prediction step k
        f(indices) = -2*C'*Qy*y_ref; % Linear term from tracking cost expansion
    end

    % GUROBI SOLVER
    model.obj = f; % Linear term
    model.rhs = [bineq; B_eq]; % Constraints RHS

    params.OutputFlag = 0; 
    result = gurobi(model,params); % Solve Gurobi

    if ~isfield(result,'x') % stop if no solution is returned
        disp(result.status);
        return;
    end

    x.opt = result.x; % Optimal decision vector [X; U]

    % Extract optimal control sequence
    u_opt = reshape(x.opt(z_end + 1:u_end), nu, n_pred);

    % Apply first control move
    u_array(:,kk) = u_opt(:,1);

    % State propagation
    if kk < n_sim
        x_array(:,kk+1) = A*x_array(:,kk) + Bu*u_array(:,kk) + Bv*v_array(:,kk); % Simulate next state
        y_array(:,kk+1) = C*x_array(:,kk+1); %  Compute zone temperatures
    end

end


%% Calculate Total Heating Demand
Sq_m = zeros(nu,1); % Allocate storage sq_m
for i = 1:nu
    identifier = sprintf('Z%04d',i);
    % get sq meter values from from zones
    Sq_m(i) = B.thermal_model_data.getValue(identifier,'area');
end

% Convert radiator heat flux [W/m^2] to total heating power [W].
% Q_tot = sum(zone area * radiator heat flux) for all zones.
Q_tot = Sq_m'*u_array;

Q_max = max(Q_tot);   % Maximum heating demand [W]
Q_min = min(Q_tot);   % Minimum heating demand [W]
Q_mean = mean(Q_tot); % Average heating demand [W]

% Total delivered heating energy [kWh].
% Power [W] * time [h] = energy [Wh].
E_total_kWh = sum(Q_tot)*Ts_hrs/1000; % total delivered heating energy [kWh]

disp(E_total_kWh);




%% Calculate buffer tank energy requirement

% Difference between actual and average heating demand [W].
% Positive values mean demand is above average, negative values below.
Q_deviation = Q_tot - Q_mean;

% Accumulate the power deviation over time to get buffer energy profile [J].
% Ts_hrs*3600 converts each timestep from hours to seconds.
E_buffer = cumsum(Q_deviation*Ts_hrs*3600);

% Required buffer capacity is the difference between highest and lowest stored energy [J].
E = max(E_buffer) - min(E_buffer);

disp('Total Energy the tank needs to store')
disp(E)


%% Calculate tank volume

% Water tank volume from stored energy: E = rho*V*cp*DeltaT
rho_water = 1000;   % Water density [kg/m^3]
cp_water = 4186;    % Specific heat capacity of water [J/(kg*K)]
Delta_T = 20;        % Usable temperature difference [K]

V = E/(rho_water*cp_water*Delta_T); % Required tank volume [m^3]

disp('Volume:');
disp(V);

%% Alternative signed buffer calculation

% This version uses the opposite sign convention.
%
% Positive value:
% Mean demand is larger than actual demand.
%
% Negative value:
% Actual demand is larger than mean demand.
%Q_deviation_signed = Q_mean - Q_tot;

% Accumulated signed buffer energy profile.
%E_buffer_signed = cumsum(Q_deviation_signed*Ts_hrs*3600);
%% Plot and save results


figure;
plot(Q_tot, 'LineWidth', 1.5);
xlabel('Timestep');
ylabel('Q_{total} [W]');
hold on;
yline(Q_mean); % reference mean
title('Total Heating Demand');
legend('Total heating demand','Mean heating demand','Location','best');
exportgraphics(gcf, fullfile(folder, 'Tot_Heat_Demand.png'));

figure;
plot(E_buffer, 'LineWidth', 1.7);
xlabel('Timestep');
ylabel('E_{buffer} [J]');
title('Buffer tank');
grid on;
exportgraphics(gcf, fullfile(folder, 'Buffer_tank.png'));

y_plot = C*x_array; % zone temperatures

figure;
plot(y_plot');
hold on;
yline(22,'r--'); % reference temperature
title('Zone temp');
xlabel('Time step');
ylabel('Temperature [°C]');
legend('Zone1','Zone2','Zone3','Zone4','Zone5');
grid on;
exportgraphics(gcf, fullfile(folder, 'ZoneTemps.png'));

figure;
plot(u_array');
title('Control input radiator 1-5');
xlabel('Time step');
ylabel('q [W/m^2]');
legend('Rad1','Rad2','Rad3','Rad4','Rad5');
grid on;
exportgraphics(gcf, fullfile(folder, 'Rad_q_w_m2.png'));

figure;
plot(Q_deviation, 'LineWidth',1.7)
xlabel('Timestep')
ylabel('Q_{deviation} [W]')
title('Buffer tank charging and discharging')
grid on
exportgraphics(gcf, fullfile(folder, 'ChargeDCharge1_buffer.png'));

Q_deviation_signed = Q_mean - Q_tot;
E_buffer_signed = cumsum(Q_deviation_signed*Ts_hrs*3600);

figure;
plot(E_buffer_signed, 'LineWidth',1.7)
xlabel('Timestep')
ylabel('E_{buffer} [J]')
title('Charge and Discharge - Buffer tank')
grid on
exportgraphics(gcf, fullfile(folder, 'ChargeDCharge2_buffer.png'));


%% Store tank sizing results

% Store important outputs in a struct for later use.
tank_size.Q_max = Q_max;
tank_size.Q_min = Q_min;
tank_size.Q_mean = Q_mean;
tank_size.E_buffer = E;
tank_size.Q_tot = Q_tot;
tank_size.E_buffer_profile = E_buffer;
tank_size.E_total_kWh = E_total_kWh;
tank_size.V = V;

% Save tank sizing data if needed.

%save('building3_data.mat','tank_size');


