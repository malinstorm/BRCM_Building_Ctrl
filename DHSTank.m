
%% DHS tank temperature update
function TankTemp_DHS_k = DHSTank(TankTemp_DHS,Q_DHS_out,b_DHS,Q_in,T_DHS_min,T_DHS_max)
%function TankTemp_DHS_k = DHSTank(TankTemp_DHS,Q_DHS_out,b_DHS,T_DHS_min,T_DHS_max)

%Q_in = 3e5;
Q_out = sum(Q_DHS_out); % Total heat supplied from DHS to buildings [W]

TankTemp_DHS_k = TankTemp_DHS + b_DHS*(Q_in - Q_out); % Update DHS tank temperature

TankTemp_DHS_k = min(T_DHS_max, TankTemp_DHS_k); % Limit to maximum tank temperature
TankTemp_DHS_k = max(T_DHS_min, TankTemp_DHS_k); % Limit to minimum tank temperature

end