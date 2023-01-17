% Simplified simulation of scooting a chair/sled forward w/o touching ground
% Author:  Mike Kokko
% Updated: 17-Jan-2023

% restart
close all; clear; clc;

% parameters
sysParams.m_a   = 12;    % [kg] chair mass
sysParams.m_b   = 75;    % [kg] person mass
sysParams.mu_s  = 0.30;  % coefficient of static friction
sysParams.mu_k  = 0.20;  % coefficient of kinetic friction
sysParams.g     = 9.81;  % [m/s^2] acceleration of gravity
sysParams.P_fwd = 4;     % [N] force pushing person CM forward (static)
sysParams.P_rev = -265;  % [N] force flinging person CM backward (dynamic) -- this should be negative!
sysParams.A_len = 0.3;   % [m] length over which person CM moves/slides
sysParams.Fn01  = (sysParams.m_a + sysParams.m_b)*sysParams.g;
sysParams.mode = 1;
sysParams.Ff = nan;

% simulation time parameters
t0 = 0;        % [s] simulation start time
tf = 6.0;      % [s] simulation end time
dt = 0.001;    % [s] timestep size

% initial conditions
x0 = zeros(4,1);
X = x0;

% data storage
time = t0;
data = x0;
mode_data = sysParams.mode;

% confirm friction assumptions
assert(abs(sysParams.P_fwd) < sysParams.mu_s*sysParams.Fn01,'Forward force too large!');
assert(abs(sysParams.P_rev) > sysParams.mu_s*sysParams.Fn01,'Reverse force too small!');

% run simulation
for t = t0:dt:(tf-dt)

    % calculate timestep for ODE solving
    odeTime = [t t+dt];

    % extract from state vector
    x_a = X(1);
    x_a_dot = X(2);
    x_b = X(3);
    x_b_dot = X(4);

    % determine whether we need to switch modes
    % and update direction of frictional force if necessary
    switch(sysParams.mode)
        case 1
            % in a static static push forward ( static friction )
            % wait for collision
            % set into state 2:
            % - link velocities w/ cons. of linear momentum
            % - kinetic friction
            if( (x_b - x_a) >= sysParams.A_len)
                X(3) = X(1) + sysParams.A_len;
                v_new = (sysParams.m_a*x_a_dot + sysParams.m_b*x_b_dot)/(sysParams.m_a + sysParams.m_b);
                X(2) = v_new;
                X(4) = v_new;
                sysParams.mode = 2;
                %                 fprintf('switching to state 2, v_new = %0.4f m/s\n',v_new);
            end
        case 2
            % masses linked/locked sliding forward subject to kinetic friction
            % wait for motion to stop
            % set into state 3:
            % - push mass B backwards strongly
            % - kinetic friction
            if( x_b_dot < 0)
                X(2) = 0;
                X(4) = 0;
                sysParams.mode = 3;
                %                 fprintf('switching to state 3\n');
            end
        case 3
            % mass a sliding on ground while mass b forced backward
            % wait for collision
            % set into state 4:
            % - link velocities with cons. of linear momentum
            % - kinetic friction acting to slow system down
            if( (x_b - x_a) <= 0 )
                X(3) = X(1);
                v_new = (sysParams.m_a*x_a_dot + sysParams.m_b*x_b_dot)/(sysParams.m_a + sysParams.m_b);
                X(2) = v_new;
                X(4) = v_new;
                sysParams.Ff = -1*sign(v_new)*sysParams.mu_k*sysParams.Fn01;
                sysParams.mode = 4;
                %                 fprintf('switching to state 4, x_a_dot was %0.4f, v_new = %0.4f m/s\n',x_a_dot,v_new);
            end
        case 4
            % masses linked/locked sliding together subject to kinetic friction
            % wait for motion to stop
            % set into state 1:
            % - push mass A forward weakly
            % - static friction
            sysParams.Ff = -1*sign(x_a_dot)*sysParams.mu_k*sysParams.Fn01;
            if( abs(x_b_dot) < 0.001 )
                X(2) = 0;
                X(4) = 0;
                sysParams.mode = 1;
                %                 fprintf('switching to state 1\n');
            end
    end

    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()

    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
    mode_data(end+1) = sysParams.mode;
end

%% animate results
figure;
set(gcf,'Position',[0680 0751 0560 0227]);
hold on;
grid on; axis equal;
ph_a = plot([nan nan],[0 0],'-','LineWidth',4,'Color',[0 0 0.8]);
ph_b = plot(nan,0,'.','MarkerSize',50,'Color',[0.8 0 0]);
xlim ([min(data(1,:)) max(data(1,:))+sysParams.A_len]);
ylim(0.1*max(xlim)*[-1 1]);
th = title(sprintf('%8.3fs',0));

for t_idx = 1:50:length(time)
    ph_a.XData = data(1,t_idx)+[0 sysParams.A_len];
    ph_b.XData = data(3,t_idx);
    th.String = sprintf('%8.3fs',time(t_idx));
    drawnow;
end

%% plot data

% speeds
fh1 = generateModeFig(time,mode_data,min([data(2,:) data(4,:)]),max([data(2,:) data(4,:)]));
figure(fh1);
set(gcf,'Position',[0323 0374 1187 0420]);
ph1(1) = plot(time,data(2,:),'-','Color',[0 0 0.8],'LineWidth',1.6);
ph1(2) = plot(time,data(4,:),'-','Color',[0.8 0 0],'LineWidth',1.6);
xlabel('\bfTime [s]');
ylabel('\bfAbsolute Speed [m/s]');
legend('Chair','Person');
legend(ph1,{'Chair','Person'},'Location','SouthWest');

% positions
fh2 = generateModeFig(time,mode_data,min([data(1,:) data(3,:)]),max([data(1,:) data(3,:)]));
figure(fh2);
set(gcf,'Position',[0323 0374 1187 0420]);
ph2(1) = plot(time,data(1,:),'-','Color',[0 0 0.8],'LineWidth',1.6);
ph2(2) = plot(time,data(3,:),'-','Color',[0.8 0 0],'LineWidth',1.6);
xlabel('\bfTime [s]');
ylabel('\bfPosition [m]');
legend(ph2,{'Back of Chair','Person CM'},'Location','NorthWest');


% propagate state
function X_dot = stateProp(t,X,sysParams)

% deconstruct state vector
% generate system matrix
x_a     = X(1);
x_a_dot = X(2);
x_b     = X(3);
x_b_dot = X(4);

X_dot = zeros(4,1);
X_dot(1) = x_a_dot;
X_dot(3) = x_b_dot;

switch(sysParams.mode)
    case 1
        % mass a static
        X_dot(2) = 0;
        % mass b being pushed forward
        X_dot(4) = (1/sysParams.m_a)*sysParams.P_fwd;
    case 2
        % masses a and b linked
        X_dot(2) = ((-sysParams.mu_k*sysParams.Fn01)/(sysParams.m_a + sysParams.m_b));
        X_dot(4) = X_dot(2);
    case 3
        % mass a sliding forward
        X_dot(2) = (1/sysParams.m_a)*(-sysParams.mu_k*sysParams.Fn01-sysParams.P_rev);
        X_dot(4) = (1/sysParams.m_b)*(sysParams.P_rev);
    case 4
        % masses a and b linked, friction should be stopping them
        X_dot(2) = (sysParams.Ff/(sysParams.m_a + sysParams.m_b));
        X_dot(4) = X_dot(2);
end
end

% generate figure with time on x axis
% and colored bands indicating when system is in each mode
function fh = generateModeFig(time,mode_data,ymin,ymax)
mode_colors = [ 0.8 0 0; 0 0.8 0; 0 0 0.8; 0.8 0.8 0];
fh = figure;
hold on; grid on;
xlim([0 time(end)]);
ylim([ymin ymax]);
this_mode = mode_data(1);
this_mode_start = 1;
for t_idx = 1:length(time)
    if (mode_data(t_idx) == this_mode) && t_idx ~= length(time)
        continue
    end
    patch('XData',time([this_mode_start t_idx t_idx this_mode_start]),...
        'YData',[ymin ymin ymax ymax],'FaceColor',mode_colors(this_mode,:),...
        'EdgeColor','none','FaceAlpha',0.1);

    % prepare for next patch
    this_mode = mode_data(t_idx);
    this_mode_start = t_idx;
end
end