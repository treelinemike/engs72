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
len = 1.0;           % [m] rod length
m_total = 1.0;       % [kg]
impulse_mag = 1.0;   % [kg*m/s]
mass_range = [2:20 25:5:400];  % number of masses

figure;
axh(1) = subplot(1,4,1);
hold on; grid on; axis equal;
ylim([-0.1 1.1]);

axh(2) = subplot(1,4,2:4);
hold on; grid on;
xlim([0 max(mass_range)]);
ylim([-0.1 1.1]);

linkaxes(axh,'y');

r0_data = [];
for num_masses = mass_range

    % distribute masses
    mass_locs = linspace(0,len,num_masses);

    rcm = sum(mass_locs)/num_masses;

    % compute vcm_plus
    vcm_plus = (impulse_mag)/(m_total);

    % compute moment of inertia at CM
    I = 0;
    for mass_idx = 1:length(mass_locs)
        r = mass_locs(mass_idx) - rcm;
        I = I + (m_total/num_masses)*r^2;
    end

    % compute angular speed
    omega = (-rcm)*impulse_mag/(I);

    % compute location of instant center (relative to impacted end of rod)
    r0 = (rcm)-vcm_plus/omega;
    r0_data(end+1,:) = [num_masses,r0];

    subplot(1,4,1);
    cla;
    plot([0 0],[0 1],'-','LineWidth',3,'Color',0.5*ones(1,3));
    plot(zeros(size(mass_locs)),mass_locs,'.','MarkerSize',30,'Color',[0 0 0.8]);
    subplot(1,4,2:4);
    plot(r0_data(:,1),r0_data(:,2),'-','LineWidth',2,'Color',[0 0 0.8]);
    drawnow;
    
    if(find(mass_range == num_masses) < 20)
        pause(0.4);
    end
end

% %%
% syms rcm r3
% eqns = [rcm == (1+r3)/(3), -3*rcm*(0.5-rcm) == - rcm^2 - (1-rcm)^2 - (rcm-r3)^2];
% S = solve(eqns,[r3,rcm])
% S.r3