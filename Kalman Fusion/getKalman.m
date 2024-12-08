function [estimated_altitude,altitude_gps,altitude_bar_interp,time_gps]=getKalman(barometer_data,gps_data)
% barometer_data = readtable(filename,Sheet="Pressure"); % Baro output
time_bar = barometer_data.Time_s_;    % Time from barometer data
pressure = barometer_data.Pressure_hPa_; % Pressure in hPa
pressure = pressure*100; % Change pressure to Pa

% Convert Pressure to Altitude
P0 = 101325; % Standard atmospheric pressure at sea level (Pa)
T0 = 288.15; % Standard temperature (K)
L = 0.0065;  % Temperature lapse rate (K/m)
R = 8.31432; % Universal gas constant (J/mol*K)
g0 = 9.80665; % Gravitational acceleration (m/s^2)
M = 0.0289644; % Molar mass of Earth's air (kg/mol)

altitude_bar = (T0/L) * (1 - (pressure / P0).^(R*L/(g0*M)));

% Load GPS data
% gps_data = readtable(filename,Sheet="Location"); % Location output
time_gps = gps_data.Time_s_;        % Time from GPS data
altitude_gps = gps_data.Height_m_; % GPS altitude in meters

R_e = 6371000; % Radius of Earth (m)
vertical_uncertainty = gps_data.VerticalAccuracy_m_; % Vertical Accuracy (m)
%vertical_uncertainty = vertical_uncertainty * pi/180 * R_e; % Convert from deg to m

% Interpolate barometer data to GPS timestamps for synchronization

%%%%%
altitude_bar=altitude_bar+altitude_gps(1)-altitude_bar(1);
%%%%%
altitude_bar_interp = interp1(time_bar, altitude_bar, time_gps, 'linear');

% Kalman Filter Initialization
dt = mean(diff(time_gps)); % Average time step based on GPS time

F = [1 dt; 0 1]; % State transition matrix
H = [1 0];       % Measurement matrix (both sensors measure altitude)

Q = [0.1 0; 0 0.01]; % Process noise covariance
R_bar = 2;           % Barometer noise covariance

x = [altitude_bar_interp(1); 0]; % Initial state (altitude and velocity)
P = eye(2); % Initial covariance

% Kalman Filter Loop
estimated_altitude = zeros(length(time_gps), 1);

for k = 1:length(time_gps)
    % Prediction Step
    x = F * x;
    P = F * P * F' + Q;
    
    % GPS noise covariance
    R_gps = vertical_uncertainty(k)^2;

    % Measurement Update (Barometer)
    z_bar = altitude_bar_interp(k); % Interpolated barometer altitude
    y_bar = z_bar - H * x;
    S_bar = H * P * H' + R_bar;
    K_bar = P * H' / S_bar;
    x = x + K_bar * y_bar;
    P = (eye(2) - K_bar * H) * P;

    % Measurement Update (GPS)
    z_gps = altitude_gps(k); % GPS altitude
    y_gps = z_gps - H * x;
    S_gps = H * P * H' + R_gps;
    K_gps = P * H' / S_gps;
    x = x + K_gps * y_gps;
    P = (eye(2) - K_gps * H) * P;

    % Store estimated altitude
    estimated_altitude(k) = x(1);
end
end
% % Plot Results
% figure;
% plot(time_gps, altitude_bar_interp, 'b-', 'DisplayName', 'Interpolated Barometer Altitude');
% hold on;
% plot(time_gps, altitude_gps, 'r*', 'DisplayName', 'GPS Altitude');
% plot(time_gps, estimated_altitude, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Estimated Altitude');
% xlabel('Time (s)');
% ylabel('Altitude (m)');
% legend;
% title('Kalman Filter: Altitude Estimation from Barometer and GPS');
% grid on;

