function app_visualizeSimulation(app, states, waves, xVec, yVec, tVec, faces, vertices, cogVec)
% 船只的状态（位置、姿态等）和波浪高度数据结合起来，生成一个动态的三维可视化效果。
% VISUALIZESIMULATION Visualize 3D simulation of a ship states on given waves.

% 清除之前的图形
cla(app.simualte_UIAxes);

% patch函数创建船只的三维模型，并设置其颜色、光照等属性。faces三角网格的面，vertices顶点坐标
p = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', [0.8 0.8 1.0], ... % 船的颜色为淡蓝色
         'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
         'AmbientStrength', 0.35, 'Parent', app.simualte_UIAxes); 
hold on;

% axis('image', app.simualte_UIAxes);
axis(app.simualte_UIAxes, 'image');
% 设置坐标轴的比例为“图像模式”，即保持 x 和 y 轴的比例一致
view(app.simualte_UIAxes, [-135 35]); % 三维视角
camlight left; % 添加左侧光源
lighting phong; % 使用Phong光照模型
material('dull'); % 设置材质为哑光

% sea = surf(app.simualte_UIAxes, xVec, yVec, waves(:, :, 1)); % surf(X,Y,Z)将矩阵 Z 中的值绘制为由 X 和 Y 定义的 x-y 平面中的网格上方的高度。
% 曲面的颜色根据 Z 指定的高度而变化。
sea = surf(app.simualte_UIAxes, xVec, yVec, waves(:, :, 1));
colormap(app.simualte_UIAxes, abyss); 
colorbar(app.simualte_UIAxes); % 添加颜色条

% 坐标轴范围，限制可视化区域
axis(app.simualte_UIAxes, [5 300 0 100 -10 35]); % IMPORTANT: change these axes for visualization of corridor 根据场景可调
xlabel(app.simualte_UIAxes, 'x [m]');
ylabel(app.simualte_UIAxes, 'y [m]');
zlabel(app.simualte_UIAxes, 'z [m]');
title(app.simualte_UIAxes, '3D Simulation Result');
grid on;
drawnow;

deltaX = 0; deltaY = 0; deltaZ = 0;

% 动态更新
for t = 2:length(tVec)-1
    sea.ZData = waves(:, :, t);

    deltaX = deltaX + states(1, t) - states(1, t-1);
    deltaY = deltaY + states(2, t) - states(2, t-1);
    deltaZ = deltaZ + states(3, t) - states(3, t-1);

     % 检查 p 是否有效
        if ~ishandle(p)
            % 如果 p 无效，重新创建它
            p = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', [0.8 0.8 1.0], ...
                     'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
                     'AmbientStrength', 0.35, 'Parent', app.simualte_UIAxes); hold on;
        end
    
    
    p.Vertices = (R(states(7, t), -states(8, t), -states(9, t))' * (vertices - cogVec(t, :))')' ...
                 + cogVec(t, :) + [deltaX -deltaY -deltaZ];
    drawnow;
end
end