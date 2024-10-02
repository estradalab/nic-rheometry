%                    Quasi-static NIC Simulator
%
%                           MAIN DRIVER
%
%   Zhiren Zhu (zhiren@umich.edu)
%   Oct. 2024
%
% =========================================================================
% USAGE
%
% Simulate the time history of pressure in needle-cavitated bubble 
% surrounded by a hyperelastic material.
%
% This script is the main driver. 
% User provides input regarding: 
%   (1) dimensional & loading info of experiment simulated, 
%   (2) choice of hyperelastic material model type and parameter range 
% The input info is then passed to a solver to obtain the corresponding
% stress and pressure history.
%
% The experiment is assumed to be quasi-static. The governing equation is
% the Rayleigh-Plesset equation with inertial terms removed. We seek an
% equilibrium between the bubble pressure, surface tension, and elastic
% stress integral.
%
% References:
%       (1) Raayai-Ardakani, S., Chen, Z., Earl, D. R., & Cohen, T. (2019). 
%           Volume-controlled cavity expansion for probing of local elastic 
%           properties in soft materials. Soft matter, 15(3), 381-392. 
%               * Cohen Group's NIC experiment 
%       (2) Yang, J., Cramer III, H. C., & Franck, C. (2020). Extracting 
%           non-linear viscoelastic material properties from 
%           violently-collapsing cavitation bubbles. Extreme Mechanics 
%           Letters, 39, 100839.
%               * Quadratic hyperelastic model is introduced here in IMR
%               setting (i.e., with inertial term and wave propagation)
%
% =========================================================================

clearvars;
clc; close all;

%% User input

% (A) Experiment setting 
R0 = 2E-4;          % (m) Equilibrium (initial) radius of bubble
Rmax = 40*R0;       % (m) Maximum radius of bubble
BA = inf;            % Ratio of outer vs. inner radius of undeformed sample

nstep = 200;         % # of time steps to simulate

% Note: the solution accuracy does not depend on step size. The choice of 
%       'nstep' is purely for result presentation.

% (B) Hyperelastic model
% (B-1) Model type
model = 'quad';

% Options:
%       'nh'       Neo-Hookean
%       'quad'     Quadratic Law

% (B-2) Parameter range:
GX     = 10.^(4);                    % Array of elastic shear modulus (Pa) 
alpX   = 10.^([-inf,-2:0.5:-1]);       % Array of strain stiffening param 'alpha' (dimensionless)

% Note: alpha is not used for Neo-Hookean. When the parameter is set to 0, the
%       quadratic law converegs to Neo-Hookean.

%% Simulate

% Loading history:
l0 = 1.0; % Initial hoop stretch
lmax = Rmax/R0; % Max hoop stretch
dl = (lmax-l0)/(nstep-1);

LAM = l0:dl:lmax; 

% Set up data storage:
nG = length(GX);
nalp = length(alpX);

PX = zeros(nG,nalp,nstep); % History of bubble pressure
SX = zeros(nG,nalp,nstep); % History of stress integral

for ii = 1:nG
    for jj = 1:nalp
        
        params = [ GX(ii), alpX(jj) ];

        [P_here, S_here] = NIC_solver(R0, LAM, model, params, BA);
        
        PX(ii,jj,:) = P_here;
        SX(ii,jj,:) = S_here;

    end
end

%% Plot

% Adjust this section according to what you wish to plot.
% Default setting compares different alpha values. Results are scaled by
% elastic modulus.

cmap = viridis(nalp + 1); % Plus 1 to avoid bright yellow

lw = 2.0; % Line Width
ms = 6.0; % Marker Size (if used)
ftsz1 = 24; % Font Size #1
ftsz2 = 18; % Font Size #2

p_inf = 101325; % Confirm that this matches the solver 

figure(100);
tcl = tiledlayout(1,2);

% Add title
if BA < inf
    ttl_str = ["$B_0/A_0 = ",BA,"$"];
else
    ttl_str = "$B_0/A_0 \to \infty$";
end
title(tcl,join(ttl_str),'FontSize',ftsz1,'Interpreter','latex')

ifix = 1;

for kk = 1:2
    nexttile;
    hold on; grid on; box on;
    pbaspect([1,1,1]);
    set(gca,'TickLabelInterpreter','Latex','FontSize',ftsz2)
end

for jj = 1:nalp
    nexttile(1);
    plot(LAM, reshape(SX(ifix,jj,:),[nstep,1])/GX(ifix), '-', ...
        'LineWidth', lw, 'MarkerSize', ms, 'Color',cmap(jj,:), ...
        'DisplayName',strcat('$\alpha=',num2str(alpX(jj),2),'$'));

    nexttile(2);
    plot(LAM, (reshape(PX(ifix,jj,:),[nstep,1]) - p_inf)/GX(ifix),'-', ...
        'LineWidth', lw, 'MarkerSize', ms, 'Color',cmap(jj,:), ...
        'DisplayName',strcat('$\alpha=',num2str(alpX(jj),2),'$'));
end

% Add labels
nexttile(1);
xl1 = xlabel('$\Lambda = R/R_0$');
set(xl1,'Interpreter','Latex','FontSize',ftsz2);
yl1 = ylabel('$S/G$');
set(yl1,'Interpreter','Latex','FontSize',ftsz2);
xlim([l0,lmax]);

nexttile(2);
xl2 = xlabel('$\Lambda = R/R_0$');
set(xl2,'Interpreter','Latex','FontSize',ftsz2);
yl2 = ylabel('$\left( p_{\rm b} - p_{\infty} \right)/G$');
set(yl2,'Interpreter','Latex','FontSize',ftsz2);
xlim([l0,lmax]);

lg = legend;
set(lg,'Interpreter','Latex','FontSize',ftsz2,'Location','Northwest', ...
    'NumColumns',2);

set(gcf, 'Position', get(0, 'Screensize'));