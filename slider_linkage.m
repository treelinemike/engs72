% restart
close all; clear; clc;

% parameters
theta = 18*pi/180;  % [rad/s]
L = 10/12;          % [ft]
g = 32.2;           % [ft/s^2]
dt = 0.001;          % [s]
ramp_len = 30/12;      % [ft]

time_vals = 0:dt:sqrt(4*L/(g*sin(theta)));

% storage
phi_vals = [];

% initialize figure
figure;
hold on; grid on; axis equal;
xlim([0 2]);
ylim([-2 2]);
plot([0 ramp_len*cos(theta)],[0 -ramp_len*sin(theta)],'-','LineWidth',1.6,'Color',[0 0 0]);
ph_l1 = plot([0 0],[nan nan],'-','LineWidth',1.6,'Color',[0 0 0.8]);
ph_l2 = plot([0 0],[nan nan],'-','LineWidth',1.6,'Color',[0 0 0.8]);

plot(0,0,'.','MarkerSize',40,'Color',[0 0 0]);

a = g*sin(theta);

drawnow;

for t = time_vals
    s = 0.5*a*t^2;
    phi = acos(s/(2*L));
    phi_vals(end+1) = phi;
    
    ph_l1.XData = [0 L*cos(phi-theta)];
    ph_l1.YData = [0 L*sin(phi-theta)];
    ph_l2.XData = [L*cos(phi-theta) s*cos(theta)];
    ph_l2.YData = [L*sin(phi-theta) -s*sin(theta)];
    
    drawnow;


end

figure;
ax(1) = subplot(2,1,1);
hold on; grid on;
plot(time_vals,phi_vals*180/pi,'LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfphi [deg]');

ax(2) = subplot(2,1,2);
hold on; grid on;
plot(time_vals,gradient(phi_vals,time_vals),'LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfphi dot [rad/s]');

linkaxes(ax,'x');

% delta_h = (20/12)*sin(theta);
% 
% A = [0.258 1.780; 0.793 5.477];
% b = [5.477 -1.780]';
% x = A\b
% 
% wab = [0 0 x(1)]';
% wbc = [0 0 x(2)]';
% 
% rba = (10/12)*[cos(theta) -sin(theta) 0]';
% rbc = (10/12)*[-cos(theta) sin(theta) 0]';
% 
% vc = sqrt(2*32.2*delta_h)*[cos(theta) -sin(theta) 0]';
% vb1 = cross(wab,rba)
% vb2 = vc + cross(wbc,rbc)
