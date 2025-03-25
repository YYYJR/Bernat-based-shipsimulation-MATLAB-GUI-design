function visualizeSimulation(states, waves, xVec, yVec, tVec, faces, vertices, cogVec)
%船只的状态（位置、姿态等）和波浪高度数据结合起来，生成一个动态的三维可视化效果。
% VISUALIZESIMULATION Visualize 3D simulation of a ship states on given waves.
figure;
%patch函数创建船只的三维模型，并设置其颜色、光照等属性。faces三角网格的面，vertices顶点坐标
p = patch('Faces', faces,'Vertices', vertices, 'FaceColor', [0.8 0.8 1.0], ...%船的颜色为淡蓝色
         'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
         'AmbientStrength', 0.35); hold on;
axis('image');%设置坐标轴的比例为“图像模式”，即保持 x 和 y 轴的比例一致
view([-135 35]);%三维视角
% camlight('headlight');%头顶光
% material('dull');%哑光材质
camlight left; % 添加左侧光源
lighting phong; % 使用Phong光照模型
material('dull'); % 设置材质为哑光

% sea = surf(xVec, yVec, waves(:, :, 1)); %surf(X,Y,Z)将矩阵 Z 中的值绘制为由 X 和 Y 定义的 x-y 平面中的网格上方的高度。
% %曲面的颜色根据 Z 指定的高度而变化。
% colorbar;  % 添加颜色条

sea = surf(xVec, yVec, waves(:, :, 1));
colormap(abyss); 
colorbar; % 添加颜色条


%坐标轴范围，限制可视化区域
axis([5 300 0 100 -10 35]); % IMPORTANT: change these axes for visualization of corridor 根据场景可调
xlabel('x [m]');
ylabel('y [m]');
pause(0.1);
deltaX = 0; deltaY = 0; deltaZ = 0;

%动态更新
for t=2:length(tVec)-1
    sea.ZData = waves(:, :, t);

    deltaX = deltaX + states(1, t) - states(1, t-1);
    deltaY = deltaY + states(2, t) - states(2, t-1);
    deltaZ = deltaZ + states(3, t) - states(3, t-1);

    p.Vertices = (R(states(7,t), -states(8,t), -states(9,t))' * (vertices-cogVec(t,:))')'...
                 + cogVec(t,:) + [deltaX -deltaY -deltaZ];
    drawnow;
end
end
