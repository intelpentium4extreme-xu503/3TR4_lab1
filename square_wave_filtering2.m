clc;
clear all;
hold off;

f0 = 4000;     %fundamental freq of input square wave
T0 = 1/f0;  %period 
tstep = 0.005*T0;
no_sample = 3*T0/tstep + 1; %no. of samples  within  3*T0
no_sample1 = T0/tstep + 1; %no. of samples  within  T0

tt = -1.5*T0:tstep:1.5*T0;

amplitude = 5; 
gp1 = amplitude * square(2 * pi * f0 * tt); 
gp_in = gp1; % Input Waveform

figure(1)
Hp1 = plot(tt,gp_in);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('input - time domain')


%% Fourier series representation of signal (Amplitude Spectrum)
      
% K=1/(2*pi);
N=100; %no. of harmonics
nvec = -N:N;
c_in = zeros(size(nvec));
for n = nvec
    m = n+N+1;
    
    if mod(n, 2) == 1
        c_in(m) = amplitude * (2 / (1i * pi * n)); % Fourier coefficients for square wave
    else
        c_in(m) = 0; % Zero coefficients for even harmonics
    end
end
f = nvec*f0; %frequency vector
figure(2)
Hp1=stem(f,abs(c_in));
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of input')


%% Fourier series representation of signal (Phase Spectrum)

figure(3)
Hp1=stem(f,angle(c_in));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
axis([-0.1e4 0.1e4 -pi pi])
title('phase spectrum of input')

%% Designing the 2nd order Butterworth filter

R= 2e3;
C=0.1e-6;
fc = 1/(2*pi*R*C);    %cutoff freq of filter

Q = 5; % Quality factor for Butterworth filter
f_gain = 3 - (1 / Q); % Filter Gain


Hf = 1 ./ (1 - (f / fc).^2 + 1i * (f / (fc * Q)));
c_out = c_in .* Hf; %Fourier coefficients of the filter output

figure(4)
stem(f,abs(c_in),'r','LineWidth',2);
hold on
stem(f,abs(c_out),'b','LineWidth',2);
hold off
axis([-8*f0 8*f0 0 max(abs(c_in))])
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of filter output and input')
Ha = gca;
set(Ha,'Fontsize',16)
legend('input','output')



%% Construct the output signal from the Cout Fourier coefficients

A = zeros(2*N+1,ceil(no_sample));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*tt);
end
gp_out = sum(A);

% Limit Vout in range  (filter_gain * Vin, -filter_gain * Vin)
V_limit = amplitude * f_gain;  % Limiting factor
gp_out = max(min(real(gp_out), V_limit), -V_limit);

figure(5)
Hp1 = plot(tt,real(gp_out),'b',tt,gp_in,'r');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('filter input and output-time domain')
set(Ha,'Fontsize',16)
legend('output','input')