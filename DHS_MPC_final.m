%function Q_DHS = DHS_MPC(T_DHS, Qdhs_out_demand, b_DHS, Qmax, T_DHS_min, T_DHS_max,lambda,kappa)
function [Q_DHS, Q_in] = DHS_MPC_final(T_DHS, Q_req, b_DHS, Qmax, Q_in_max, T_DHS_min, T_DHS_max, lambda, kappa)

n_buildings = length(Q_req); % Number of connected buildings
n_pred = 24; % Prediction horizon

n_u = n_buildings + 1; % Number of control inputs: Q_out for each building + Q_in to DHS tank
n_dec = n_u*n_pred; % Total decision variables over prediction horizon

indices_Q_out = 1:n_buildings; % Indices for heat supplied to buildings ([q1,q2,q3,q4])
indices_Q_in  = n_buildings + 1; % Index for heat production into DHS tank (qin, place 5 ctrl vec)

%% Cost

rho_DHS = 1e-10; % Small regularization term for numerical stability
H = 2*rho_DHS*eye(n_dec); % Quadratic cost matrix

lambda = lambda(:); % Dual price signal from coupling algorithm

prod_cost = 0.195; %0.1;%0.02; % [kr/kWh]
alpha = 1e3; % weighting factor for cost lambda - not used

f = zeros(n_dec,1); % Linear cost vector

for k = 1:n_pred
    j_start = (k-1)*n_u;

    Q_out_indices = j_start + indices_Q_out; % Heat supplied to buildings
    Q_in_indices  = j_start + indices_Q_in;  % Heat produced into DHS tank
    % putting it in linear cost
    f(Q_out_indices) = -lambda(:)*(0.25/1000); % Revenue from supplying heat
    f(Q_in_indices)  = prod_cost*(0.25/1000); % Production cost for DHS heat
end




%% Ineq Constraints

% 1) DHS supply capacity: sum Q_out <= Qmax
% The total heat supplied to all buildings cannot exceed Qmax
Aineq_capacity = zeros(n_pred, n_dec); % allocation storage
bineq_capacity = Qmax*ones(n_pred,1);

for k = 1:n_pred % Loop over prediction horizon

    j_start = (k-1)*n_u; % Start index for prediction step k
    Aineq_capacity(k, j_start + indices_Q_out) = 1; % Set coefficients for Q_out terms in the capacity constraint(used later by Gurobi in A*x <= b)

end




% 2) Minimum DHS tank temperature: T_DHS(k) >= T_DHS_min
Aineq_temp_min = zeros(n_pred, n_dec); % Allocate constraint matrix
bineq_temp_min = zeros(n_pred,1); % RHS

for k = 1:n_pred % Constraint at prediction step k

    for j = 1:k % Accumulate all previous heat flows affecting tank temperature

        j_start = (j-1)*n_u;
        % Heat supplied to buildings lowers tank temperature
        Aineq_temp_min(k, j_start + indices_Q_out) = Aineq_temp_min(k, j_start + indices_Q_out) + b_DHS; 
         % Heat production increases tank temperature
        Aineq_temp_min(k, j_start + indices_Q_in) = Aineq_temp_min(k, j_start + indices_Q_in) - b_DHS;

    end

    bineq_temp_min(k) = T_DHS - T_DHS_min;% RHS of minimum temperature constraint

end

% 3) Maximum DHS tank temperature: T_DHS(k) <= T_DHS_max
Aineq_temp_max = zeros(n_pred, n_dec); % Allocate constraint matrix
bineq_temp_max = zeros(n_pred,1); % RHS

for k = 1:n_pred % Constraint at prediction step k

    for j = 1:k % Include all previous heat inputs and outputs up to prediction step k

        j_start = (j-1)*n_u;
        % Heat supplied to buildings lowers tank temperature
        Aineq_temp_max(k, j_start + indices_Q_out) = Aineq_temp_max(k, j_start + indices_Q_out) - b_DHS; 
        % Heat production increases tank temperature

        Aineq_temp_max(k, j_start + indices_Q_in) = Aineq_temp_max(k, j_start + indices_Q_in) + b_DHS; 
    end

    bineq_temp_max(k) = T_DHS_max - T_DHS; % RHS of max tank temp constraint

end

% goes into gurobi as ineq constraints
Aineq = [Aineq_capacity; Aineq_temp_min; Aineq_temp_max];
bineq = [bineq_capacity; bineq_temp_min; bineq_temp_max];


%% Bounds

lb = zeros(n_dec,1); % Lower bounds: Q_out >= 0, Q_in >= 0
ub = inf(n_dec,1);   % Initialize upper bounds

for k = 1:n_pred

    j_start = (k-1)*n_u;

    ub(j_start + indices_Q_out) = Qmax; % Upper limit on heat supplied to buildings: Q_out <= Qmax
    ub(j_start + indices_Q_in)  = Q_in_max; % Upper limit on DHS heat production: Q_in <= Q_in_max

end

%% Gurobi

model.Q = sparse(H); % Hessian matrix
model.obj = f; % Linear objective function
model.A = sparse(Aineq); % Inequality constraint matrix
model.rhs = bineq; % Right-hand side of constraints
model.sense = repmat('<', size(Aineq,1),1); % Constraint type: A*x <= b
model.lb = lb; % Lower bounds
model.ub = ub; % Upper bounds

params.OutputFlag = 0; 

result = gurobi(model,params); % Solve optimization problem

Q_opt = result.x; % Optimal decision vector

%% First step of opt
Q_DHS = Q_opt(indices_Q_out); % supply to buildings
Q_in  = Q_opt(indices_Q_in); % supply to the DHS

end % end function

% Parts of the code documentation, alternative lambda-based solution approaches were developed with assistance from ChatGPT
% and subsequently reviewed, modified and validated by the author.

