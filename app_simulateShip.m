function [states, faces, vertices, cogVec] = app_simulateShip(app,wavesFile, shipStruct, isPlot, isVisual)
% SIMULATESHIP Ship on sea simulation. Given a waves file and a ship,
% simulates 12 states of the ship through time. The states are: 12个状态如下
%   - [x,     y,    z]:    position of the ship’s center of gravity (cog)
%                          described in the earth-fixed coordinate
%                          system;船舶质心在xyz轴的位置
%   - [v_u,   v_v,  v_w]:  velocity of the ship’s cog described in the body
%                          -fixed coordinate system;速度（相对船舶质心 船体固定坐标系）
%   - [phi,   th,   psi]:  Euler angles describing rotation of the body-
%                          fixed coordinate system in relation to the earth
%                          -fixed coordinate system;横摇、纵摇、艏摇角
%   - [w_phi, w_th, w_psi]: rotational velocity in the body-fixed
%                           coordinate system.上面对应的角速度
% The forces considered are the hidrostatic (buoyancy) force by the water, 
% the ship's weight, the coriolis effect and the propellers force. Ignores:
% added mass, wind forces, dynamic water forces, etc. To calculate the 
% inertia, the ship was assumed to be a solid cuboid. 
%
% Inputs:
%   - wavesFile:  file containing waves and its properties;海浪文件
%   - shipStruct: struct containing STL file and other ship properties;船结构体
%   - isPlot:     boolean, if true: plots states through time;是否绘制状态信息
%   - isVisual:   boolean, if true: show visualization of simulation in 3D.是否可视化
% Outputs:
%   - states:    all 12 states simulated through the time vector defined in
%                the wave file;12个状态信息,包含船只状态的矩阵，每一行表示一个时间步的船只位置和姿态（例如，states(1, :) 表示 x 位置，states(7, :) 表示横摇角等）
%   - faces:     ship's faces, only returned to be able to show
%                visualization outside function;船的三角面信息，用于可视化
%   - vertices:  ship's vertices, only returned to be able to show
%                visualization outside function;船舶的顶点信息
%   - cogVec:    center of gravity vector, only returned to be able to show
%                visualization outside function.船舶质心的轨迹
%
% Constants:
%   - Ax: cross-sectional area of ship along ship's x-axis [m^2];沿x轴方向的横截面积
%   - Ay: cross-sectional area of ship along ship's y-axis [m^2];
%   - Az: cross-sectional area of ship along ship's z-axis [m^2];
%   - Ix: ship's momentum of inertia along x axis;船舶绕x轴的转动惯量
%   - Iy: ship's momentum of inertia along y axis;
%   - Iz: ship's momentum of inertia along z axis;
%   - Ki: integrational constant for the PI-regulator;PI调节器的积分常数
%   - Kp: proportional constant for the PI-regulator;PI调节器的比例常数
%   - B:  damping constant for ship's rotational velocity;船舶旋转速度的阻尼系数
%   - Cd: damping constant for ship's velocity;船舶速度的阻尼系数。
%   - M:  ship's mass [kg];船的质量
%   - ro: water density [kg/m^3];水的密度
%   - g:  gravitational acceleration [m/s^2];重力加速度
%   - l:  length of ship (along x-axis) [m];船的长
%   - w:  width of ship (along y-axis) [m];宽
%   - h:  height of ship (along z-axis) [m].高

tic;

fprintf('\nStarted shipOnSea simulation!\n');

% ------------------------------- Define constants ------------------------
% ----- Non-tunable constants

% Ship dimensions inheritant to HMS Norfolk 从stl3D文件里读取数据，暂时不知道怎么改成用户输入
len    = shipStruct.len; 
width  = shipStruct.width; 
height = shipStruct.height;
M      = shipStruct.M;
% Cross sectional area of the submerged hull 横截面积
Au     = width * 0.9 * height * 0.33;
Av     = height * 0.33 * len;
Aw     = len * width * 0.9;

% Offset to ship's center of gravity, set by trial and error 
cogOffset = shipStruct.cogOffset;

% Ship's moment of inertia 转动惯量
I  = M/12 * diag([width^2 + height^2, len^2 + height^2, len^2 + width^2]);
Iu = I(1, 1);
Iv = I(2, 2); 
Iw = I(3, 3);

% Physical constants
ro = 997; 
g  = 9.81;   

% ----- Tunable constants PI调节器的比例参数
Kp_force  = 9.5;
Ki_force  = 0.82;
Kp_torque = 3.5;
Ki_torque = 0.5;

C_du = 2.5; %速度阻尼系数
C_dv = 2.5;
C_dw = 2.5;

B_phi = 1e10 / 20; %旋转阻尼
B_th  = 1e10;
B_psi = 1e10;

% ------------------------------- Load waves ------------------------------
% disp('1) Loading waves...');
waves = load(wavesFile).wavesStruct.waves;
beta  = load(wavesFile).wavesStruct.beta;
xVec  = load(wavesFile).wavesStruct.xVec;
yVec  = load(wavesFile).wavesStruct.yVec;
tVec  = load(wavesFile).wavesStruct.tVec;
Ts    = load(wavesFile).wavesStruct.Ts;
% displayName = load(wavesFile).wavesStruct.displayName;
% fprintf(['Waves loaded: \n', displayName]);

% ------------------------------- Load ship's STL file --------------------
% disp('2) Loading ship and computing normal vectors to triangles...');
[faces, vertices, ~] = stlreadOwn(shipStruct.file);

vertices = vertices + shipStruct.verticesPos; 
% Compute center of gravity (cog) from vertices (geometric center) 计算质心
cog = [(max(vertices(:, 1)) + min(vertices(:, 1))) / 2, ...
       (max(vertices(:, 2)) + min(vertices(:, 2))) / 2, ...
       (max(vertices(:, 3)) + min(vertices(:, 3))) / 2];
cog = cog + cogOffset; %质心偏移量

[facePoints, faceAreas, normals] = calculatePointsAreasNormals(vertices); %计算每个三角面的中心点、面积和法向量

disp('Boat loaded!');

% ------------------------------- Define time update model and simulate ---
% disp('3) Simulating states through time...');

% ----- Define update matrices A and B (C is arbitrary) and discretize.
% A是状态更新矩阵，B是输入矩阵，输入对状态的影响
%    x y z      v_u       v_v     v_w         phi th psi w_phi   w_th w_psi 
A = [0 0 0       1         0       0          0  0   0    0       0      0;
     0 0 0       0         1       0          0  0   0    0       0      0;
     0 0 0       0         0       1          0  0   0    0       0      0;
     0 0 0 -C_du*ro*Au/M   0       0          0  0   0    0       0      0;
     0 0 0       0 -C_dv*ro*Av/M   0          0  0   0    0       0      0;
     0 0 0       0         0 -C_dw*ro*Aw/M    0  0   0    0       0      0;
     0 0 0       0         0       0          0  0   0    1       0      0;
     0 0 0       0         0       0          0  0   0    0       1      0;
     0 0 0       0         0       0          0  0   0    0       0      1;
     0 0 0       0         0       0          0  0   0 -B_phi/Iu  0      0;
     0 0 0       0         0       0          0  0   0    0     -B_th/Iv 0;
     0 0 0       0         0       0          0  0   0    0       0  -B_psi/Iw]; 
 
%    F_x/M   F_y/M   F_z/M   Tau_x/Ix   Tau_y/Iy   Tau_z/Iz
B = [  0       0       0        0          0          0;  % x dot
       0       0       0        0          0          0;  % y dot
       0       0       0        0          0          0;  % z dot
       1       0       0        0          0          0;  % v_u dot
       0       1       0        0          0          0;  % v_v dot
       0       0       1        0          0          0;  % v_w fot
       0       0       0        0          0          0;  % phi dot
       0       0       0        0          0          0;  % th dot
       0       0       0        0          0          0;  % psi dot
       0       0       0        1          0          0;  % w_phi dot
       0       0       0        0          1          0;  % w_th dot
       0       0       0        0          0          1]; % w_psi dot
C = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0];
sys = c2d(ss(A, B, C, []), Ts);

% ----- Define reference states and initialize other variables 定义参考状态和初始化变量
refSpeedU      = getVelMetPerTs(shipStruct.refSpeedU, Ts); % m/Ts
refYaw         = shipStruct.refYaw;
extraUForce    = 0; 
extraYawTorque = 0;
cogVec         = zeros(length(tVec), 3);

% ----- Initialize all states and set 1st state
states = zeros(12, length(tVec));
%                    [x       y       z     ], [v_u v_v v_w phi th psi w_phi w_th w_psi]'
states(:, 1) = cat(2,[cog(1)  -cog(2)  -cog(3)],shipStruct.x0)'; 

% ----- Rotate facePoints & normals according to initial states  根据初始姿态调整船舶的面中心点和法向量
phi = states(7, 1); th = states(8, 1); psi = states(9, 1);
facePoints = (R(phi, -th, -psi)' * (facePoints - cog)')' + cog; 
normals    = (R(phi, -th, -psi)' * normals')'; 

% ----- State update for all time steps 更新状态变量
for tIdx=1:length(tVec)-1
    cogVec(tIdx, :) = cog;
    % ------- Compute sum of all forces F_net & sum of all torques Tau_net
    % 计算力和力矩
    F_net = zeros(3, 1);
    Tau_net = zeros(3, 1);
    for nIdx=1:length(normals)
        % Get index of waves that is closest to the evaluated normal
        xIdx = round(facePoints(nIdx, 1)); 
        yIdx = round(facePoints(nIdx, 2)); 
        facePointHeight = facePoints(nIdx, 3);
        % 限制索引范围
        xIdx = max(1, min(size(waves, 2), xIdx));
        yIdx = max(1, min(size(waves, 1), yIdx));
        % Add netto component if wave is above the considered normal
        if waves(yIdx, xIdx, tIdx) > facePointHeight   %如果海浪高于当前面的高度，计算浮力和力矩    
            % Force
            h = waves(yIdx, xIdx, tIdx) - facePointHeight;                                    
            volume = h * faceAreas(nIdx);
            F_buoy = (ro * g * volume) * normals(nIdx, :)';
            F_net = F_net + F_buoy;      
            % Torque
            lever = (facePoints(nIdx, :) - cog)';
            phi = states(7, tIdx); th = states(8, tIdx); psi = states(9, tIdx);
            Tau = cross(R(phi, -th, -psi) * lever, R(phi, -th, -psi) * F_buoy);
            Tau_net = Tau_net + Tau;
        end
    end
    
    % ------- Set up input
    phi = states(7, tIdx); th = states(8, tIdx); psi = states(9, tIdx);
    v_u = states(4, tIdx);
    % 使用PI调节器计算额外的推力和转矩
    extraUForce = Ki_force * extraUForce + pRegulator(Kp_force, refSpeedU, v_u);
    forcesInLocalCoord = R(phi, -th, -psi) * [F_net(1)/M; -F_net(2)/M; -F_net(3)/M + g] ... 
                         + [extraUForce; 0; 0];
    
    extraYawTorque = Ki_torque * extraYawTorque + pRegulator(Kp_torque, refYaw, psi);
    torqueInLocalCoord = [Tau_net(1)/Iu; -Tau_net(2)/Iv; -Tau_net(3)/Iw] ...
                         + [0; 0; extraYawTorque];
    
    inputs = [forcesInLocalCoord' torqueInLocalCoord']';
    
    % ------- Time-update 状态更新
    v_u= states(4, tIdx); v_v= states(5, tIdx); v_w= states(6, tIdx);
    w_phi = states(10, tIdx); w_th = states(11, tIdx); w_psi = states(12, tIdx);
    velInGlobalCoord = R(phi, th, psi)' * [v_u; v_v; v_w];
    rotVelDerivative = T(phi, th) * [w_phi; w_th; w_psi];
    
    coriolisV = cross([w_phi, w_th, w_psi]', [v_u, v_v, v_w]');
    
    % Make time-update as below due to non-linearities in the A matrix
    states(:, tIdx+1) = [states(1,tIdx) + sys.A(1,4) * velInGlobalCoord(1);
                         states(2,tIdx) + sys.A(2,5) * velInGlobalCoord(2);
                         states(3,tIdx) + sys.A(3,6) * velInGlobalCoord(3);
                         sys.A(4,4) * states(4,tIdx) - Ts * coriolisV(1);
                         sys.A(5,5) * states(5,tIdx) - Ts * coriolisV(2);
                         sys.A(6,6) * states(6,tIdx) - Ts * coriolisV(3);
                         states(7,tIdx) + sys.A(7,10) * rotVelDerivative(1);
                         states(8,tIdx) + sys.A(8,11) * rotVelDerivative(2);
                         states(9,tIdx) + sys.A(9,12) * rotVelDerivative(3);
                         sys.A(10,10) * states(10,tIdx);
                         sys.A(11,11) * states(11,tIdx);
                         sys.A(12,12) * states(12,tIdx)] + sys.B * inputs;
    
    % ------- Update ship hull 更新几何信息
    deltaX = states(1, tIdx+1) - states(1, tIdx);
    deltaY = states(2, tIdx+1) - states(2, tIdx);
    deltaZ = states(3, tIdx+1) - states(3, tIdx);
    
    deltaPhi = states(7, tIdx+1) - states(7, tIdx);
    deltaTh = states(8, tIdx+1) - states(8, tIdx);
    deltaPsi = states(9, tIdx+1) - states(9, tIdx);
    
    % Rotate facepoints and normals with R'
    facePoints = (R(deltaPhi, -deltaTh, -deltaPsi)' * (facePoints - cog)')' ...
                 + cog + [deltaX -deltaY -deltaZ];
    normals = (R(deltaPhi, -deltaTh, -deltaPsi)' * normals')';
    cog = cog + [deltaX -deltaY -deltaZ];
end
disp('Simulation done!');
toc;

% ------------------------------- Plot states through time ----------------
if isPlot
    app_plotShipStates(app,states, tVec, Ts, beta);%绘制状态变量随时间的变化
%     plotHelipadStates(states, tVec, Ts, beta, shipStruct.helipadPos);
end
% ------------------------------- Visualize simulation in 3D --------------
if isVisual
    app_visualizeSimulation(app,states, waves, xVec, yVec, tVec, faces, vertices, cogVec);
end

%% Help functions 前面运用到的辅助函数定义
function [facePoints, faceAreas, normals] = calculatePointsAreasNormals(V)
    % Returns facepoints which are the center of the triangles, the
    % triangles' areas and the normals to the triangles.
    j = 1;
    for i=1:3:length(V)-2
        facePoints(j, :) = [mean(V(i:i+2,1)), mean(V(i:i+2, 2)), mean(V(i:i+2, 3))];
        faceAreas(j) = areaOfFace(V(i:i+2,1), V(i:i+2, 2), V(i:i+2, 3)); 
        p0 = V(i, :)'; 
        p1 = V(i+1, :)'; 
        p2 = V(i+2, :)';
        n = cross(p0-p1, p0-p2)' ./ norm(cross(p0-p1, p0-p2)'); 
        normals(j, :) = -n; % The normal of the face
        j = j + 1;
    end
    facePoints(isnan(facePoints)) = 0;
    normals(isnan(normals)) = 0;
end
function T = T(phi, th)
% Gustafsson, "Statistical Sensor Fusion" 3rd edition
% Eq. (13.9), p. 349
    T = [1 sin(phi)*tan(th) cos(phi)*tan(th);
         0       cos(phi)         -sin(phi);
         0 sin(phi)/cos(th) cos(phi)/cos(th)];
end
function regulation = pRegulator(Kp, refValue, currentValue)
% Simple P-control of speed.
    regulation = Kp * (refValue - currentValue);    
end
function velMeterPerTs = getVelMetPerTs(vel, Ts)
% Gets a velocity in m/s and returns a velocity in m/Ts
    velMeterPerTs = vel * Ts;
end

end