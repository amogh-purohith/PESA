% ----- Parameters -----
c = 3e8;                % speed of light (m/s)
f = 5e9;                % carrier frequency (Hz)
lambda = c / f;         % wavelength (m)
M = 4;                  % elements in x-direction
N = 4;                  % elements in y-direction
S_e = 1;                % element pattern (isotropic)

% Observation grid
theta_deg = -90:2:90;
phi_deg   = -180:2:180;
theta = deg2rad(theta_deg);
phi   = deg2rad(phi_deg);

% Element spacing
buf1 = 0.010;           % 1 cm buffer
d = lambda/2 - buf1;
dx = d; dy = d;
k = 2*pi / lambda;

% Precompute observation grid (fixed)
[Theta, Phi] = meshgrid(theta, phi);
u = sin(Theta) .* cos(Phi);
v = sin(Theta) .* sin(Phi);
X = sin(Theta) .* cos(Phi);
Y = sin(Theta) .* sin(Phi);
Z = cos(Theta);
mask = Z >= 0;                         % keep only upper hemisphere
X(~mask) = NaN; Y(~mask) = NaN; Z(~mask) = NaN;

% Sweep angles (serpentine)
theta_s_vals = -90:10:0;               % outer loop: elevation, step 10°
phi_s_vals   = -90:5:90;               % inner loop: azimuth, step 5°
theta_s_rad = deg2rad(theta_s_vals);
phi_s_rad   = deg2rad(phi_s_vals);

% Set up figure with two subplots side by side
fig = figure('Position', [100 100 1400 600]);

% ---- Left subplot: Normalized Power (linear) ----
ax_left = subplot(1,2,1);
hold on; grid on;
view([-90,45,90]);
xlabel('x = sin\theta cos\phi');
ylabel('y = sin\theta sin\phi');
zlabel('z = cos\theta');
title('Normalized Power Pattern (Linear)');
axis equal;
clim([0 1]);
colormap(ax_left, turbo);
lighting gouraud;

% ---- Right subplot: Gain Pattern (dB) ----
ax_right = subplot(1,2,2);
hold on; grid on;
view([-90,45,90]);
xlabel('x = sin\theta cos\phi');
ylabel('y = sin\theta sin\phi');
zlabel('z = cos\theta');
title('Gain Pattern (dB)');
axis equal;
clim([-40 0]);               % dB range: 0 dB (peak) down to -40 dB
colormap(ax_right, turbo);
lighting gouraud;

% Initial calculation for first steering angle (to create surface objects)
u_s0 = sin(theta_s_rad(1)) * cos(phi_s_rad(1));
v_s0 = sin(theta_s_rad(1)) * sin(phi_s_rad(1));
si_x0 = k * dx * (u - u_s0) + 1e-12;
si_y0 = k * dy * (v - v_s0) + 1e-12;
S_x0 = (exp(1j * si_x0 * (M - 1)) .* sin(0.5 * M * si_x0)) ./ sin(0.5 * si_x0);
S_y0 = (exp(1j * si_y0 * (N - 1)) .* sin(0.5 * N * si_y0)) ./ sin(0.5 * si_y0);
S_total0 = S_e .* S_x0 .* S_y0;
p0 = abs(S_total0).^2;
epsVal = 1e-12;
p_dB0 = 10*log10(p0 + epsVal) - 10*log10(max(p0(:)) + epsVal);
p_n0 = p0 / max(p0(:));

% Mask invalid points
p_n0(~mask) = NaN;
p_dB0(~mask) = NaN;

% Create surface objects in both subplots
h_left  = surf(ax_left,  X, Y, Z, p_n0,  'EdgeColor', 'none');
h_right = surf(ax_right, X, Y, Z, p_dB0, 'EdgeColor', 'none');

% Text object for dynamic angles (shared, placed above both subplots)
angleText = annotation('textbox', [0.42 0.92 0.16 0.05], ...
    'String', '', 'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'BackgroundColor', 'none', ...
    'EdgeColor', 'none');

% Serpentine sweep (no video recording)
for i = 1:length(theta_s_rad)
    theta_s = theta_s_rad(i);
    theta_s_deg = theta_s_vals(i);

    % Determine phi order: forward for odd i, backward for even i
    if mod(i, 2) == 1
        phi_indices = 1:length(phi_s_rad);          % ascending
    else
        phi_indices = length(phi_s_rad):-1:1;       % descending
    end

    for j_idx = 1:length(phi_indices)
        j = phi_indices(j_idx);
        phi_s = phi_s_rad(j);
        phi_s_deg = phi_s_vals(j);

        % Update the dynamic angle text
        set(angleText, 'String', sprintf('\\theta_s = %d°, \\phi_s = %d°', theta_s_deg, phi_s_deg));

        % Direction cosines of the source
        u_s = sin(theta_s) * cos(phi_s);
        v_s = sin(theta_s) * sin(phi_s);

        % Array factor for current steering direction
        si_x = k * dx * (u - u_s) + 1e-12;
        si_y = k * dy * (v - v_s) + 1e-12;
        S_x = (exp(1j * si_x * (M - 1)) .* sin(0.5 * M * si_x)) ./ sin(0.5 * si_x);
        S_y = (exp(1j * si_y * (N - 1)) .* sin(0.5 * N * si_y)) ./ sin(0.5 * si_y);
        S_total = S_e .* S_x .* S_y;

        % Power and normalized quantities
        p = abs(S_total).^2;
        p_n = p / max(p(:));               % linear normalized power
        p_dB = 10*log10(p + epsVal) - 10*log10(max(p(:)) + epsVal);  % dB gain

        % Apply hemisphere mask
        p_n(~mask) = NaN;
        p_dB(~mask) = NaN;

        % Update both plots
        set(h_left,  'CData', p_n);
        set(h_right, 'CData', p_dB);
        drawnow;
    end
end
