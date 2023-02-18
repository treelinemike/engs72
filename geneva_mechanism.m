% restart
close all; clear; clc;

% parameters
a = 90;          % [mm]
b = 50;          % [mm]
phi_dot = 2.240; % [rad/s]
phi_0 = -pi/2;   % [rad]
phi_f = pi/2;    % [rad]
dt = 0.001;      % [s]

% construct time
t = 0:dt:((phi_f-phi_0)/phi_dot);
phi_vals = t*phi_dot + phi_0;

% storage
theta_vals = [];
theta_dot_analytical_mak = [];
theta_dot_analytical_team = [];


%%
% setup figure
figure;
hold on; grid on; axis equal;
xlim([0 90]);
ylim(1.25*b*[-1 1]);
ph_pb = plot(nan(1,2),nan(1,2),'-','LineWidth',3,'Color',[ 0.8 0.4 0 ]);
ph_oa = plot(nan(1,2),nan(1,2),'.-','MarkerSize',40,'LineWidth',3,'Color',[ 0 0.8 0 ]);
plot(0,0,'.','MarkerSize',40,'Color',[0 0 0]);
plot(a,0,'.','MarkerSize',40,'Color',[0 0 0]);

for phi = phi_vals

    c = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
    theta = asin((b/c)*sin(phi));
    theta_vals(end+1) = theta;

    theta_dot_analytical_team(end+1) = ( phi_dot*cos(phi + theta) / (cos(theta)-cos(phi + theta)) );
%     theta_dot_analytical_mak(end+1) = (2*(b^2)*sin(phi)*cos(phi)*phi_dot - (sin(theta)^2)*2*a*b*sin(phi)*phi_dot )/(2*sin(theta)*cos(theta)*(a^2+b^2-2*a*b*cos(phi)));
    theta_dot_analytical_mak(end+1) = (-b*phi_dot*cos(pi-theta-phi))/(a*cos(theta)+b*cos(pi-theta-phi));


    ph_oa.XData = [0 b*cos(phi)];
    ph_oa.YData = [0 b*sin(phi)];
    ph_pb.XData = [a a-c*cos(theta)];
    ph_pb.YData = [0 c*sin(theta)];
%     drawnow;
    
end

%%
figure;
subplot(2,1,1);
hold on; grid on;
plot(t,phi_vals);
plot(t,theta_vals);

subplot(2,1,2);
hold on; grid on;
plot(t,gradient(phi_vals,t));
plot(t,gradient(theta_vals,t));

figure;
ax(1) = subplot(3,1,1);
hold on; grid on;
plot(phi_vals*180/pi,theta_vals*180/pi,'-','LineWidth',1.6);
xlabel('\bfphi [deg]');
ylabel('\bfTheta [deg]');

ax(2) = subplot(3,1,2:3);
hold on; grid on;
plot(phi_vals*180/pi,gradient(phi_vals,t),'-','LineWidth',1.6);
plot(phi_vals*180/pi,gradient(theta_vals,t),'-','LineWidth',1.6);
plot(phi_vals*180/pi,theta_dot_analytical_team,'--','LineWidth',1.6);
plot(phi_vals*180/pi,theta_dot_analytical_mak,'--','LineWidth',1.6);
xlabel('\bfphi [deg]');
ylabel('\bfAngular Speed [rad/s]');
legend('phi dot','theta dot','team analytical','mike analytical','Location','NorthWest');
ylim([-1 6]);

linkaxes(ax,'x');
