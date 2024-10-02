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
% This script is the solver.
% =========================================================================

function [PB, SI] = NIC_solver(R0, LAM, model, params, BA)

% Input:
%   R0: Initial radius of bubble (m)
%   LAM: history of hoop stretch
%   params: array containing viscoelastic parameters
%   model:  type of hyperelastic model
%   BA: Ratio of outer vs. inner radius of undeformed sample

% Output:
%   PB: history of bubble pressure
%   SI: history of stress integral

% Read out material parameters
G = params(1);        % (Pa) Elastic shear modulus
alp = params(2);      % (dimensionless) strain stiffening parameter

% Additional parameters:
gam = 0.072;        % (N/m) Surface tension
p_inf = 101325;     % (Pa) Atmospheric pressure

% Surface tension term:
ST = -2*gam./(R0*LAM);

% Get outer wall stretch:
BLAM = ((LAM.^3 -1)*(BA^(-3)) + 1).^(1/3);

% Now, calculate history
switch model
    case 'nh'
        SI = -0.5*G*(BLAM.^(-4) + 4./BLAM - LAM.^(-4) - 4./LAM);
        
    case 'quad'
        chunk1 = (1/8)*BLAM.^(-8) + (1/5)*BLAM.^(-5) + BLAM.^(-2) - 2*BLAM;
        chunk2 = (1/4)*BLAM.^(-4) + BLAM.^(-1);
        chunk3 = (1/8)*LAM.^(-8) + (1/5)*LAM.^(-5) + LAM.^(-2) - 2*LAM;
        chunk4 = (1/4)*LAM.^(-4) + LAM.^(-1);
        SI = 2*G*((3*alp-1)*(chunk2 - chunk4) - alp*(chunk1 - chunk3));
    otherwise
        error('Invalid model type. Please reselect.');
end

% Bubble pressure:
PB = p_inf - SI - ST;

end