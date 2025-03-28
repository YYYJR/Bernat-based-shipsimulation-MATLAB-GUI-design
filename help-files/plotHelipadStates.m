function plotHelipadStates(states, tVec, Ts, beta, helipadPos)
% PLOTSHIPSTATES Plots the 12 states of the helipad through time vector 
% tVec.
%这段代码定义了一个函数 `plotHelipadStates`，用于绘制船（helipad）的状态随时间的变化。
%它基于输入的船舶状态、时间向量、采样时间间隔、波浪方向和船甲板的位置，计算并绘制甲板的全局位置、速度、姿态角和角速度。


figure;    
r = states(1:3, :) - helipadPos';

for i = 1:length(states)
   globalHelipadPos(:,i) = R(states(7,i),states(8,i),states(9,i))*r(:,i);

   localVel = cross([states(10,i),states(11,i),states(12,i)]', r(:,i))+states(4:6,i);
   globalVel(:,i) = R(states(7,i),states(8,i),states(9,i))*localVel;
end

% Global position states x, y, z
subplot(4, 3, 1); plot(tVec * Ts, globalHelipadPos(1, :)); title('Global Helipad position X [m]'); xlabel('Time [s]');
subplot(4, 3, 2); plot(tVec * Ts, globalHelipadPos(2, :)); title('Global Helipad position Y [m]'); xlabel('Time [s]');
subplot(4, 3, 3); plot(tVec * Ts, -globalHelipadPos(3, :)); title('Global Helipad position Z (UP)[m]'); xlabel('Time [s]');

% Local velocity states u, v, w
subplot(4, 3, 4); plot(tVec * Ts, globalVel(1, :) / Ts); title('Global velocity U [m/s]'); xlabel('Time [s]');
subplot(4, 3, 5); plot(tVec * Ts, globalVel(2, :) / Ts); title('Global velocity V [m/s]'); xlabel('Time [s]');
subplot(4, 3, 6); plot(tVec * Ts, -globalVel(3, :) / Ts); title('Global velocity W (UP)[m/s]'); xlabel('Time [s]');

% Local angles phi, th, psi
subplot(4, 3, 7); plot(tVec * Ts, rad2deg(states(7, :))); title('Roll \phi [deg]'); xlabel('Time [s]');
subplot(4, 3, 8); plot(tVec * Ts, rad2deg(states(8, :))); title('Pitch \theta [deg]'); xlabel('Time [s]');
subplot(4, 3, 9); plot(tVec * Ts, rad2deg(states(9, :))); title('Yaw \psi [deg]'); xlabel('Time [s]');

% Local rotational velocities w_phi, w_th, w_psi
subplot(4, 3, 10); plot(tVec * Ts, rad2deg(states(10, :)) / Ts); title('Roll velocity w_\phi [deg/s]'); xlabel('Time [s]');
subplot(4, 3, 11); plot(tVec * Ts, rad2deg(states(11, :)) / Ts); title('Pitch velocity w_\theta [deg/s]'); xlabel('Time [s]');
subplot(4, 3, 12); plot(tVec * Ts, rad2deg(states(12, :)) / Ts); title('Yaw velocity w_\psi [deg/s]'); xlabel('Time [s]');

sgtitle(['Helipad states. Long crested waves from ', num2str(rad2deg(beta)), ' degrees']);
end