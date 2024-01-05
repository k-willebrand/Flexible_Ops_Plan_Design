%% create T_Temp for deterministic path
function T_Temp = deterministicTempMatrix(N);

T_Temp_single = zeros(N);
trans = eye(N-1);

T_Temp_single(2:N,1:N-1) = trans;
T_Temp_single(end,end) = 1;

T_Temp = repmat(T_Temp_single, [1 1 5]);