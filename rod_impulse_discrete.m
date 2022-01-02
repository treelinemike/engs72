% discrete approximation of delivering impulse to end of uniform rod
% and computing the position of the instant center just after impulse is
% applied
%
% when using two masses (one at each end of rod) instant center is at the
% far end of the rod
%
% when using infinite number of masses, instant center is at 2/3 length of
% rod away from the impacted end (as we know from hw problem)

% restart
close all; clear; clc;

% parameters
num_masses = 10000;  % number of masses to distribute
len = 1.0;           % [m] rod length
m_total = 1.0;       % [kg]
impulse_mag = 1.0;   % [kg*m/s]

% compute mass locations
mass_locs = linspace(0,len,num_masses);

% compute vcm_plus
vcm_plus = (impulse_mag)/(m_total);

% compute moment of inertia at CM
I = 0;
for mass_idx = 1:length(mass_locs)
    r = (len/2) - mass_locs(mass_idx);
    I = I + (m_total/num_masses)*r^2;
end

% compute angular speed
omega = (-len/2)*impulse_mag/(I);

% compute location of instant center (relative to impacted end of rod)
r0 = (len/2)-vcm_plus/omega
