function app_plotShipStates(app, states, tVec, Ts, beta)
% PLOTSHIPSTATES Plots the 12 states of the ship through time vector tVec.
% 用于绘制船舶的 12 个状态变量随时间的变化。这些状态变量包括船舶的全局位置、局部速度、姿态角和角速度。
% 函数的输入参数包括船舶的状态矩阵、时间向量、采样时间间隔和波浪方向角。

% 时间向量转换为秒
timeInSeconds = tVec * Ts;

% 绘制全局位置
plot(app.GCG_X, timeInSeconds, states(1, :));
title(app.GCG_X, 'Global CoG position X [m]');
xlabel(app.GCG_X, 'Time [s]');

plot(app.GCG_Y, timeInSeconds, states(2, :));
title(app.GCG_Y, 'Global CoG position Y [m]');
xlabel(app.GCG_Y, 'Time [s]');

plot(app.GCG_Z, timeInSeconds, -states(3, :));
title(app.GCG_Z, 'Global CoG position Z (UP) [m]');
xlabel(app.GCG_Z, 'Time [s]');

% 绘制局部速度
plot(app.vU, timeInSeconds, states(4, :) / Ts);
title(app.vU, 'Local velocity U [m/s]');
xlabel(app.vU, 'Time [s]');

plot(app.vV, timeInSeconds, states(5, :) / Ts);
title(app.vV, 'Local velocity V [m/s]');
xlabel(app.vV, 'Time [s]');

plot(app.vW, timeInSeconds, -states(6, :) / Ts);
title(app.vW, 'Local velocity W (UP) [m/s]');
xlabel(app.vW, 'Time [s]');

% 绘制姿态角
plot(app.R, timeInSeconds, rad2deg(states(7, :)));
title(app.R, 'Roll \phi [deg]');
xlabel(app.R, 'Time [s]');

plot(app.P, timeInSeconds, rad2deg(states(8, :)));
title(app.P, 'Pitch \theta [deg]');
xlabel(app.P, 'Time [s]');

plot(app.Y, timeInSeconds, rad2deg(states(9, :)));
title(app.Y, 'Yaw \psi [deg]');
xlabel(app.Y, 'Time [s]');

% 绘制角速度
plot(app.wR, timeInSeconds, rad2deg(states(10, :)) / Ts);
title(app.wR, 'Roll velocity w_\phi [deg/s]');
xlabel(app.wR, 'Time [s]');

plot(app.wP, timeInSeconds, rad2deg(states(11, :)) / Ts);
title(app.wP, 'Pitch velocity w_\theta [deg/s]');
xlabel(app.wP, 'Time [s]');

plot(app.wY, timeInSeconds, rad2deg(states(12, :)) / Ts);
title(app.wY, 'Yaw velocity w_\psi [deg/s]');
xlabel(app.wY, 'Time [s]');

% % 设置总标题
% sgtitle(app.boatstate_part, ['Long crested from ', num2str(rad2deg(beta)), ' degrees']);
end