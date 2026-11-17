clear; clc; close all;

%% Physical Constants & Magnet Parameters
rng('shuffle'); 
mu0 = 4*pi*1e-7;
gamma = 1.760859e11;
kB = 1.380649e-23;
alpha = 0.01;
Ms = 1.1e6;
T = 300;
diameter = 10e-9;
thickness = 1e-9;
radius = diameter/2;
V = pi*radius^2*thickness;

% Calculate Anisotropy for Low Barrier Magnet
Ku = (2.5 * kB * T) / V; 
Hk = 2 * Ku / (mu0 * Ms);
Delta_kT = (Ku * V) / (kB * T);
disp(['Anisotropy field = ', num2str(Hk), ' A/m']);
disp(['Energy barrier = ', num2str(Delta_kT), ' kBT']);

% No external field, zero spin current
h_ext = [0; 0; 0];
hs = 0;
p = [0; 0; 1]; 

%% Time setup & Monte Carlo Parameters
dt = 1e-13;
tEnd = 15e-9; 
dtau = dt * gamma * mu0 * Hk / (1 + alpha^2);
tauEnd = tEnd * gamma * mu0 * Hk / (1 + alpha^2);
tau_arr = 0:dtau:tauEnd;
N = length(tau_arr);
time = tau_arr * (1 + alpha^2) / (gamma * mu0 * Hk);

num_runs = 1000;           % Total simulation runs requested by professor
sample_interval = 5;       % Sample every 5th time instant

% Discard initial 3 ns as burn-in thermalization
t_burn_in = 3e-9; 
burn_in_idx = find(time >= t_burn_in, 1, 'first');

D = Ms / Hk; 

% Initialize pool to aggregate mz data from all runs
all_mz = [];
m_last = zeros(N, 3);      % To store the trajectory of the final run

%% Monte Carlo Integration Loop (1000 Runs)
fprintf('Starting %d simulation runs with multi-instant sampling...\n', num_runs);
tic;
for r = 1:num_runs
    theta0 = 0.01;        % Near mz = +1
    phi0 = 2 * pi * rand;
    
    m = zeros(N, 3);
    m(1,:) = [sin(theta0)*cos(phi0), sin(theta0)*sin(phi0), cos(theta0)];
    
    % Thermal Field Amplitude (Normalized for dtau)
    sigma = sqrt((2 * alpha * kB * T) / (mu0 * Ms * V * Hk * (1 + alpha^2) * dtau));
    
    % Heun's Method Integration
    for n = 1:N-1
        m0 = m(n,:)';
        h_th = sigma * randn(3, 1);
        
        k1 = llg_sun_rhs(m0, alpha, h_ext, h_th, hs, p, D);
        mp = m0 + dtau * k1;
        mp = mp / norm(mp);
        
        k2 = llg_sun_rhs(mp, alpha, h_ext, h_th, hs, p, D);
        mnew = m0 + 0.5 * dtau * (k1 + k2);
        mnew = mnew / norm(mnew);
        
        m(n+1,:) = mnew';
    end
    
    % Capture multiple time instants after the burn-in period
    run_mz = m(burn_in_idx:sample_interval:end, 3);
    all_mz = [all_mz; run_mz(:)];
    
    % Save the last run's full trajectory for time-series plotting
    if r == num_runs
        m_last = m;
    end
    
    % Print progress update every 100 runs
    if mod(r, 100) == 0
        fprintf('Completed %d / %d runs...\n', r, num_runs);
    end
end
toc;

%% Plotting Results
plot_results(time, m_last, all_mz);

%% LLG RHS (Easy Axis Z, Hard Axis X)
function dm = llg_sun_rhs(m, alpha, h_ext, h_th, hs, p, D)
    m = m/norm(m);
    
    % Uniaxial Anisotropy along Z (Easy Axis)
    h_ani = [0; 0; m(3)];
    
    % Optional demagnetization / shape effect components
    Nx = 1; Ny = 0; Nz = 0; 
    h_demag = -D * [Nx*m(1); Ny*m(2); Nz*m(3)];
    
    h_eff = h_ani + h_demag + h_ext + h_th;
    
    t_prec = -cross(m, h_eff);
    t_damp = -alpha * cross(m, cross(m, h_eff));
    t_stt  = -hs * cross(m, cross(m, p));
    
    dm = t_prec + t_damp + t_stt;
end

%% Plotting Function (Original Style)
function plot_results(time, m_last, all_mz)
    mz_last = m_last(:,3);
    
    figure('Color', 'w');
    plot(time*1e9, mz_last, 'k', 'LineWidth', 1.2);
    xlabel('Time (ns)', 'FontSize', 12);
    ylabel('m_z (Easy Axis)', 'FontSize', 12);
    title('Telegraph Noise along Easy Axis (Z)');
    grid on; box on; ylim([-1.1, 1.1]);
    
    % PDF of m_z showing symmetric U-shaped distribution around +1 and -1
    figure('Color', 'w');
    histogram(all_mz, 100, 'Normalization', 'pdf');
    xlabel('m_z', 'FontSize', 12);
    ylabel('PDF', 'FontSize', 12);
    title('Symmetric PDF for Easy-Axis Z Magnet');
    grid on; box on;
end