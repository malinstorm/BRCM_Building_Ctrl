%% Initial DHS water tank setup
Ts_hrs = 0.25; % Sampling time [h]
n_days = 1; % Simulation length [days]
n_sim = n_days*24*4; % Number of simulation steps

%% DHS tank parameters
V_DHS = 11; % DHS tank volume [m^3]
rho = 1000; % Water density [kg/m^3]
cp = 4186; % Specific heat capacity of water [J/(kg·K)]
Ts = Ts_hrs*3600; % Sampling time [s]

b_DHS = Ts/(rho*cp*V_DHS); % Tank temperature conversion factor [°C/W]

%% DHS tank state
TankTemp_DHS = zeros(1,n_sim+1); % DHS tank temperature trajectory
TankTemp_DHS(1) = 70; % Initial DHS tank temperature [°C]