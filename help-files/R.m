function R = R(phi, th, psi)
% R Rotation matrix using the xyz rule.
% Gustafsson, "Statistical Sensor Fusion" 3rd edition
% Eq. (13.7), p. 348
R = [cos(th)*cos(psi)       cos(th)*sin(psi)       -sin(th);
     -cos(phi)*sin(psi)+sin(phi)*sin(th)*cos(psi) cos(phi)*cos(psi)+sin(phi)*sin(th)*sin(psi) sin(phi)*cos(th);
     sin(phi)*sin(psi)+cos(phi)*sin(th)*cos(psi) -sin(phi)*cos(psi)+cos(phi)*sin(th)*sin(psi) cos(phi)*cos(th)];
end

%`R` 函数根据输入的欧拉角生成旋转矩阵，适用于将向量从局部坐标系转换到全局坐标系。
%转坐标系