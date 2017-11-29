% addSpecie - Incoherent Scatter Radar Spectrum Simulator for MATLAB
% -------------------------------------------------------------------------
% Function to Add Species to the Incoherent Scatter Radar simulator for
% multiple Ion species
% -------------------------------------------------------------------------

% To determine the Atomic Weight of species use http://www.science.co.il/PTelements.asp
% And see reference paper: Dalipi and Syla (2013)
% Theoretical Analysis of ISR Spectra Behavior for Different Maxwellian Ionospheric Plasma Conditions

function [configuration, index]=addSpecie_ISR_simulator(configuration,Name_s)
% Define generic constants
q_e = 1.60217656535e-19;     % Electron charge [C]
amu = 1.66053904020e-27;     % Atomic Mass Unit [kg]

% Define Species Mass
mass_O  = 15.9994*amu;
mass_N  = 14.0067*amu;
mass_H  = 1.0079*amu;
mass_He = 4.0026*amu;
mass_molecularIons = 30.5*amu;

%Determine the specie Mass and Charge depending on Species name
if strcmpi(Name_s,'O+')
    charge = q_e;
    mass   = mass_O;
elseif strcmpi(Name_s,'O2+')
    charge = q_e;
    mass   = 2*mass_O;
elseif strcmpi(Name_s,'N+')
    charge = q_e;
    mass   = mass_N;
elseif strcmpi(Name_s,'N2+')
    charge = q_e;
    mass   = 2*mass_N;
elseif strcmpi(Name_s,'H+')
    charge = q_e;
    mass   = mass_H;
elseif strcmpi(Name_s,'He+')
    charge = q_e;
    mass   = mass_He;
elseif strcmpi(Name_s,'NO+')
    charge = q_e;
    mass   = mass_N+mass_O;
elseif strcmpi(Name_s,'molecularIons')
    charge = q_e;
    mass   = mass_molecularIons;
else
    error('ERROR: Unknown Specie name: %s\n',Name_s)
end

% Add a new specie to the struct
configuration.species.Name_s{end+1} = Name_s;      % The species name
configuration.species.q_s(end+1)    = charge;      % The charge of species [Coulomb]
configuration.species.m_s(end+1)    = mass;        % Mass of specie [AMU]

% Determine index
index = length(configuration.species.Name_s);
end