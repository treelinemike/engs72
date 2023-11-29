% Demonstration illustrating that coefficient of restitution
% scales the component of relative velocity normal to the plane of contact
%
% This is a planar motion scenario but could be extended to 3D in a
% straightforward manner. For one example see
% Computer Animation : Algorithms and Techniques by Rick Parent (but note
% that one of the relative velocity equations has a typo)
%
% Author: Mike Kokko
% Modified: 29-Nov-2023

% restart
close all; clear; clc;

%% Parameters for the planar motion scenario
% (two particle rigid baton falling vertically w/o initial rotation, particle A to contact flat surface)
m = 0.057*2;              % [kg]
e = 0.75;                 % coefficient of restitution
R = 1;                    % [m] overall length of baton
phi = 0*pi/180;           % [rad] INITIAL INCLINATION OF BATON FROM HORIZONTAL
vcm_1y = -sqrt(2*9.81*1); % [m/s] pre-impact vertical speed
vcm_1 = [0; vcm_1y; 0];
r_cm_to_A = (R/2)*[cos(phi); -sin(phi); 0];

%% Compute vertical linear impulse I required for collision to be conservative
% then use that value of I to compute a final state,
% double check that we did conserve energy,
% and look at resulting velocity of contact point (particle A)
% vertical component should be equal and opposite of incoming velocity
% b/c coeff. of restitution scales normal component of relative velocity

% kinetics
I_cons = (-2*m*vcm_1(2))/((cos(phi)^2)+1);
I = [0; I_cons; 0];
vcm_2 = vcm_1 + (I/m);  % = (I/m + vcm_1)
angimp = cross(r_cm_to_A,I);    
omega = angimp/(0.25*m*R^2);  % = 2*cos(phi)*I(2)/(m*R)
 
% kinematics
va_1 = vcm_1;
va_2 = vcm_2 + cross(omega,r_cm_to_A);
vb_2 = vcm_2 + cross(omega,-r_cm_to_A);

% are we still rigid?
u = unitvec(r_cm_to_A);
assert(abs(dot(va_2,u) - dot(vb_2,u)) < 10*eps,'RIGID CONSTRAINT FAILED!');

% is the additional impulse required (from other mass) directed along baton?
addl_imp = (m/2)*(va_2-va_1)-I;
if(norm(addl_imp) > 10*eps)
    ang_of_addl_impulse_on_a_deg = atand(addl_imp(2)/addl_imp(1));  % should be -phi... impulse delivered from particle B on A along axis of baton
    assert(abs((ang_of_addl_impulse_on_a_deg)*pi/180 - (-1*phi)) < 10*eps,'Additional impulse not directed along rod!');
end

% did we conserve energy? and compute energy correctly?
E1 = 0.5*m*norm(vcm_1)^2;
E2 = 0.5*m*norm(vcm_2)^2 + 0.125*m*(R^2)*(omega(3))^2; 
E2_sys = 0.5*(m/2)*norm(va_2)^2 + 0.5*(m/2)*norm(vb_2)^2;
assert(abs(E1 - E2) < 10*eps,'Energy not conserved!');
assert(abs(E2 - E2_sys) < 10*eps,'Energy computation in state 2 mismatch!');
va_2_conservative = va_2

%% Now try to repeat this using the specificed coefficient of restitution
% will only match the above result when e = 1

% solve linear system for vcm_2y, theta_dot, and va_2y
% the system comes from: 
% eq1: angular impulse-momentum about contact point
% eq2: kinematics linking va_2 to vcm_2
% eq3: impact (coefficient of restitution)
A = [   cos(phi),     -0.5*R,         0; ...
        1,            0.5*R*cos(phi), -1; ...
        0,            0,               1];
b = [vcm_1y*cos(phi); 0; -e*vcm_1y];
x = A\b;
vcm_2y = x(1);
theta_dot = x(2);
va_2y = x(3);

% additional kinematics
va_2x = 0.5*R*theta_dot*sin(phi);
vcm_2 = [0; vcm_2y; 0];
va_1 = vcm_1;
va_2 = [va_2x; va_2y; 0]
omega = [0; 0; theta_dot];
vb_2 = vcm_2 + cross(omega,-r_cm_to_A);

% are we still rigid?
u = unitvec(r_cm_to_A);
assert(abs(dot(va_2,u) - dot(vb_2,u)) < 10*eps,'RIGID CONSTRAINT FAILED!');

% is the additional impulse required (from other mass) directed along baton?
I = m*(vcm_2-vcm_1);
addl_imp = (m/2)*(va_2-va_1)-I; % impulse from other mass on particle A
if(norm(addl_imp) > 10*eps)
    ang_of_addl_impulse_on_a_deg = atand(addl_imp(2)/addl_imp(1));  % should be -phi... impulse delivered from particle B on A along axis of baton
    assert(abs((ang_of_addl_impulse_on_a_deg)*pi/180 - (-1*phi)) < 10*eps,'Additional impulse not directed along rod!');
end