function waves = app_simulateWaves(app,seaState, xVec, yVec, beta, tVec, U, lambda, muVec, dmu)
%Bretschneider 谱模型生成海浪高度数据
% SIMULATEWAVES(seaState, xVec, yVec, beta, tVec, U , lambda, muVec, dmu) 
% Takes in sea state and plots a wave height. Uses the Bretschneider spectrum.海浪谱
%
% Inputs:
%   - seaState: integer in interval [1, 9].
%   - xVec:     1xN vector. E.g. xVec = linspace(-50, 50, 100).
%   - yVec:     1xN vector. E.g. yVec = linspace(-50, 50, 100).
%   - beta:     direction of main wave in rad.波的传播方向
%   - tVec:     time vector. E.g.: tVec = 0:0.1:50.
%   - U:        speed of the ship [m/s]
%   - muVec:    angle directions vector. E.g.: -pi/2:dmu:pi/2. If muVec=[],
%               then the waves generated will be unidirectional (long-
%               crested).
%   - dmu:      direction interval taken in muVec. Used if muVec ~= [].
%   - lambda:   waveLength. If a sea with infinite depth is assumed, set
%               lambda to [].
%
% Ouput:
%   - waves:    if ~is3d, then size(waves) = (1, length(tVec)), where waves
%               is equal to the wave height in some point in the sea. If  
%               is3d, then size(waves) = (length(xVec), length(yVec), 
%               length(tVec)), where waves will then represent the wave
%               height of all points (x, y) over a 2D grid defined by xVec
%               and yVec over an interval of time tVec.
tic; 
disp('Creating waves...');
g = 9.81;

% Get significant wave height
seaState = app.sea_state.Value;
Hs = app_getSignificantWaveHeight(app,seaState);
% disp(['Significant wave height: ', num2str(Hs)]);


% Create Bretschneider spectrum given Hs and plot spectrum
dw = 0.1;
wVec = (dw/2:dw:3)';

% %Bretschneider模型
% A = 8.1 * 1e-3 * g^2; % constant, Equation (8.54) in Fossen
% B = 3.11 / (Hs^2);    % Equation (8.55) in Fossen
% specType = 1; % Bretschneider (@ Fossen pg 203) 
% S = wavespec(specType, [A, B], wVec, 1);
% S(1) = 0; % the first element is NaN for some reason

% JONSWAP (Vwind10,Fetch)
Fetch=150000; %150km[100,200km]
Vwind10 = app.vwind10.Value;% 设置为用户输入
specType =6;
S = wavespec(specType, [Vwind10, Fetch], wVec, 0);

disp('Wave spectrum:');
disp(S);
waves = zeros(length(yVec), length(xVec), length(tVec));

% Get the set of frequencies, directions and phases that will be used to
% generate waves for all points in the grid:
waveFrequencies = zeros(1, length(wVec));
for k=1:length(waveFrequencies)        
    waveFrequencies(k) = wVec(k) - dw/2 + dw * rand;
     disp(['Random frequency for index ', num2str(k), ': ', num2str(waveFrequencies(k))]);
end

if nargin > 7 % Short-crested wave
    waveDirections = zeros(1, length(muVec));
    for i=1:length(waveDirections)
        waveDirections(i) = muVec(i) - dmu/2 + dmu * rand;            
    end
    sizeWaveDirections = length(waveDirections);
else
    sizeWaveDirections = 1;
end

wavePhases = zeros(sizeWaveDirections, length(waveFrequencies));
for k=1:length(waveFrequencies) 
    for i=1:sizeWaveDirections
        wavePhases(k, i) = 2 * pi * rand;
    end 
end

%模拟波的高度，也需要考虑船的速度，因为是相当于在船上感受波的高度，不同船速和位置对波的高度相位都有影响。
% Get the wave heights for all coordinates (x,y) for all times in tVec
for yIdx = 1:length(yVec)
    y = yVec(yIdx);
    for xIdx=1:length(xVec)
        x = xVec(xIdx);

        % For each coordinate (x,y), sum the contribution of all the waves
        % to get the final wave amplitude at that point = sumOfWaves.
        sumOfWaves = 0;

        for k=1:length(waveFrequencies)
            w_k = waveFrequencies(k);
            
            % Long-crested
            if nargin < 8
                if nargin > 6
                    % Lambda is passed as an argument
                    coeff = 2 * pi / lambda; 
                else
                    % Infinite depth sea assumed
                    coeff = w_k ^ 2 / g;
                end
                e_k = wavePhases(k);
                amp = sqrt(2 * S(k) * dw);
                wave = amp * cos(coeff * ((x+U*tVec) * cos(-beta) ...
                               + y * sin(-beta))  ...
                                          - w_k * tVec + e_k);
                sumOfWaves = sumOfWaves + wave;
                
            % Short-crested
            else              
                for i=1:length(waveDirections)
                    e_ik = wavePhases(1, k);
                    mu_i = waveDirections(i);

                    amp = sqrt(2 * S(k) * spread(mu_i) * dw * dmu);
                    wave = amp * cos(w_k^2/g * ((x+U*tVec) * cos(mu_i - beta) ...
                                   + y * sin(mu_i - beta))  ...
                                              - w_k * tVec + e_ik);
                    sumOfWaves = sumOfWaves + wave;                
                end 
            end
        end
        waves(yIdx, xIdx, :) = sumOfWaves;
    end
end   
disp('Done creating waves!');
 % 在生成波浪数据后，添加以下调试代码
    disp('Wave data at different time steps:');
    disp(waves(:, :, 1));  % 第一个时间步的波浪数据
    disp(waves(:, :, end));  % 最后一个时间步的波浪数据
toc;


% ----------------------- Help functions ----------------------------------

function spreadFunction=spread(mu)
    if (mu >= -pi/2 && mu <= pi/2)
        spreadFunction = (2 / pi) * cos(mu)^2;
    else
        spreadFunction = 0;
    end
end
end
