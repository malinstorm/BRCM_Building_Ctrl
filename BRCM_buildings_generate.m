%% GENERATING BUILDINGS FROM BRCM - EACH BUILDING MUST LOAD UNIQUE THERMAL & EHFM DATA FROM EXCEL FILES

clear; clc; close all; % remove variables from WS, clear CMD window, closes all open firgures


%% Paths
% Add all folders inside ./fcns to the MATLAB path.
% fcns contains functions such as getOmega, getPsi,...
addpath(genpath('./fcns'));

% Get root path of the BRCM toolbox/project.
BRCMRootPath = getBRCMRootPath();

% Global debug level used by BRCM
global g_debugLvl
g_debugLvl = 0;

%% Building
demoBuilding = 'demoBuilding'; % name of building model

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



%% Draw building

% 3-D plot of Building
%B.drawBuilding;

% It is possible to control the labeling
% B.drawBuilding('NoBELabels');
% B.drawBuilding('NoLabels');
% B.drawBuilding('NoZoneLabels');

% 2-D plot of Building
%B.drawBuilding('Floorplan');


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
% Disturbance indices for ambient temperature and solar radiation
% These indices are used later to fill the disturbance vector v_array
indices_Tamb_v = getIdIndex('v_Tamb',B.building_model.identifiers.v);
indices_SolE_v = getIdIndex('v_solGlobFac_E',B.building_model.identifiers.v);
indices_SolN_v = getIdIndex('v_solGlobFac_N',B.building_model.identifiers.v);
indices_SolS_v = getIdIndex('v_solGlobFac_S',B.building_model.identifiers.v);

% Store weather disturbance indices for later use
weather_idx.Tamb = indices_Tamb_v;
weather_idx.SolE = indices_SolE_v;
weather_idx.SolS = indices_SolS_v;
weather_idx.SolN = indices_SolN_v;

%% Internal Gains 
% Get disturbance indices for internal gains in each zone
% Only the indices are stored here; the actual internal gain profiles are added later
indices_IG_zone1 = getIdIndex('v_IG_ZoneOne',B.building_model.identifiers.v);
indices_IG_zone2 = getIdIndex('v_IG_ZoneTwo',B.building_model.identifiers.v);
indices_IG_zone3 = getIdIndex('v_IG_ZoneThree',B.building_model.identifiers.v);
indices_IG_zone4 = getIdIndex('v_IG_ZoneFour',B.building_model.identifiers.v);
indices_IG_zone5 = getIdIndex('v_IG_ZoneFive',B.building_model.identifiers.v);

ig_idx.zone1 = indices_IG_zone1;
ig_idx.zone2 = indices_IG_zone2;
ig_idx.zone3 = indices_IG_zone3;
ig_idx.zone4 = indices_IG_zone4;
ig_idx.zone5 = indices_IG_zone5;


%% Radiator bounds from BRCM
% HEAT FLUX - min/max W/m^2 control input
% The allowed radiator heat flux in each zone (5).
constraintsParameters.Rad.Q_rad_Radiator1_min = 0;
constraintsParameters.Rad.Q_rad_Radiator1_max = 150; % 150

constraintsParameters.Rad.Q_rad_Radiator2_min = 0;
constraintsParameters.Rad.Q_rad_Radiator2_max = 150;

constraintsParameters.Rad.Q_rad_Radiator3_min = 0;
constraintsParameters.Rad.Q_rad_Radiator3_max = 150;

constraintsParameters.Rad.Q_rad_Radiator4_min = 0;
constraintsParameters.Rad.Q_rad_Radiator4_max = 150;

constraintsParameters.Rad.Q_rad_Radiator5_min = 0;
constraintsParameters.Rad.Q_rad_Radiator5_max = 150;


%% linear constraint matrices (Fx = state constRaint + Fu input constraint <= g)
[Fx,Fu,~,g] = B.building_model.getConstraintsMatrices(constraintsParameters);


%% Model dimensions
nz_zones = 5; % number of zones in building
nx = size(A,1); % N thermal states N x N states 
nu = size(Bu,2); % Number of ctrl inputs --> 5 ctrl variables
ny = nz_zones; % Number of measured zone temepratures
nv = size(Bv, 2); % N disturbance inputs --> 10 dist variables



%% Zone areas
% Get floor area for each zone. Used later to convert heat flux [W/m^2] to power [W].
Sq_m_zones = zeros(nu,1); % allocate storage for square meter per zone
for i = 1:nu
    identifier = sprintf('Z%04d',i); % column for zone values
    Sq_m_zones(i) = B.thermal_model_data.getValue(identifier,'area'); 
end

identifiers = B.building_model.identifiers;

% Save building model data and indices for later MPC simulations
save('building3_data.mat','B','A','Bu','Bv','Fx','Fu','g','Sq_m_zones','identifiers','weather_idx','ig_idx','nx','nu','nv','ny','nz_zones','Ts_hrs','-append');




