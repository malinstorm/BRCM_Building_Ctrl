# 2026-malin-storm-andersen
Master Thesis of Malin Storm-Andersen



Dependencies and Tools Used:

BRCM-framework
	- DOWNLOAD AND INSTALL

https://control.ee.ethz.ch/software/BRCM-Toolbox.html
https://ethz.ch/content/dam/ethz/special-interest/itet/ifa/ifa-dam/BRCM_toolbox/BRCM_Documentation-v0.95.pdf

Matlab 2022
Gurobi Solver

FCNS
C:\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03




TO GENERATE BUILDINGS FROM SCRATCH, RUN THE FILES IN ORDER 1-3 AND REPEAT FOR THE NUMBER OF BUILDINGS:
THE EXCEL SHEETS FOR THE DIFFERENT BUILDINGS WILL HAVE TO PUT IN THE FOLDER EHFM and ThermalModel in

C:\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\BuildingData\DemoBuilding

The files are provided in zip DEMOBUILDING respective folders here on GitHub.


TO DO A SIMULATION WITH ALL 4 BUILDINGS INSTANTLY, JUST JUMP TO STEP 4.
 
1.

To calculate a sufficient tank volume for the buffer tanks in the buildings, 

run BRCM_MPC_Qmax.m

with the followin parameters:

%% MPC parameters
n_days = 1;
n_sim = n_days*24*4; % n_days*24 hrs, 4 times/hour
n_pred = 24; % prediction horizon 
nz_zones = 5; % number of zones in building

Be sure to make a folder to save the tank volume and associated values:

folder = 'C:\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\Buildings\Building_n\n_day';

Save .mat - file to:
save('NAME_OF_FILE.mat','tank_size');

2.

To generate necessary values for the Whole building, run 

BRCM_buildings_generate.m 

save('NAME_OF_FILE.mat','B','A','Bu','Bv','Fx','Fu','g','Sq_m_zones','identifiers','weather_idx','ig_idx','nx','nu','nv','ny','nz_zones','Ts_hrs','-append');

3.

To do a basic tuning of the MPC with the DHS and Electricity Grid Connected, run 

BRCM_MPC_tuning.m

%% MPC parameters
n_days = 1 - 30;
n_sim = n_days*24*4; % n_days*24 hrs, 4 times/hour
n_pred = 24; % prediction horizon

folder = 'C:\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\MAKE_FOLDER\BuildingN\n_day';

%% Load saved building data
load('NAME_OF_FILE.mat')

%% Save tuning and relevant data for later


save('NAME_OF_FILE.mat','mpc_tuning','mpc_results','v_array','-append');



4. The necessary .mat files are provided so that it is possible to run all buildings instantly.
'building1_data.mat','building2_data.mat','building3_data.mat','building4_data.mat' 

Make sure 

building_data = {'building1_data.mat','building2_data.mat','building3_data.mat','building4_data.mat'};

is in the file DHS_GRID_4building.m 

To save plots -->
folder = '\MATLAB\BRCM\BRCMToolbox_v1.03\brcm\v1.03\all\BRCMToolbox_v1.03\MAKE_FOLDER\BuildingN';

Then run 

DHS_GRID_4building.m

5. Lambda implementation

DHS_Building_TankCoupling.m

DHSTank.m

DHS_Water_Tank.m

DHS_MPC_final.m




6. stresstest - turn off grid

Comment in/out parts in: DHS_Building_TankCoupling.m

Dependencies:
DHSTank.m
DHS_Water_Tank.m
DHS_MPC_final.m (tune variable prod_costs (higher, heat is more expensive))
 

7.

Implement price ceiling like suggested in NVE-report

Comment in/out parts in: DHS_Building_TankCoupling.m

Dependencies:
DHSTank.m
DHS_Water_Tank.m
DHS_MPC_final.m (tune variable prod_costs (higher, heat is more expensive))


