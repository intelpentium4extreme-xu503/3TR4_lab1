clc; clear; hold off;
%initial info construction
f0 = 4000; T0 = 1/f0; tstep = 0.005*T0;
tt = -1.5*T0:tstep:1.5*T0;
amplitude = 5;

%sqwave generator
gp_in = amplitude * square(2 * pi * f0 * tt);  
figure(1); plot(tt, gp_in, 'LineWidth', 2);
title('Input - Time Domain'); set(gca, 'FontSize', 16);

%impulse generator
N = 100; nvec = -N:N; c_in = zeros(size(nvec));
c_in(mod(nvec,2) == 1) = amplitude * (2 ./ (1i * pi * nvec(mod(nvec, 2) == 1)));
f = nvec * f0;

% Magnitude and Phase Spectrum
figure(2); stem(f, abs(c_in), 'LineWidth', 2); axis([-8*f0 8*f0 0 max(abs(c_in))]); title('Magnitude Spectrum');
figure(3); stem(f, angle(c_in), 'LineWidth', 2); axis([-1e4 1e4 -pi pi]); title('Phase Spectrum');

%impulse response
R = 2e3; C = 0.1e-6; fc = 1/(2*pi*R*C); Q = 5; f_gain = 3 - (1 / Q);
Hf = 1 ./ (1 - (f / fc).^2 + 1i * (f / (fc * Q))); c_out = c_in .* Hf;

%compairsion response
figure(4); stem(f, abs(c_in), 'r', 'LineWidth', 2); hold on; stem(f, abs(c_out), 'b', 'LineWidth', 2);
axis([-8*f0 8*f0 0 max(abs(c_in))]); legend('Input', 'Output'); title('Magnitude Spectrum of Filtered Input');


% Time-domain graph
gp_out = sum(c_out .* exp(1i * 2 * pi * nvec' * f0 * tt));
V_limit = amplitude * f_gain; gp_out = max(min(real(gp_out), V_limit), -V_limit);
figure(5); plot(tt, icsreal(gp_out), 'b', tt, gp_in, 'r', 'LineWidth', 2);
legend('Output', 'Input'); title('Filter Input and Output - Time Domain'); set(gca, 'FontSize', 16);