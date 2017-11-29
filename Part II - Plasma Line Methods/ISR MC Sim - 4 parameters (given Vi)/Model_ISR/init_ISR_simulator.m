% init_ISR_simulator 
%
%       Incoherent Scatter Radar Spectrum Simulator for MATLAB
% -------------------------------------------------------------------------
% Function to initialize the Incoherent Scatter Radar simulator for 
% multiple Ion species 
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 25 January 2016
% -------------------------------------------------------------------------



function [configuration] = init_ISR_simulator(radarFrequency, minFreq, maxFreq,numPoints,Name_s)
    
    % Define generic constants
    q_e = 1.60217656535e-19;        % Electron charge [C]
    m_e = 9.1093821545e-31;         % Electron mass [kg]
    c   = 2.99792458e8;             % Speed of light [m/s]

    % Initialize configuration struct
    configuration = {};

    % Generate struct for configurations
    configuration.radarFrequency = radarFrequency;          % Radar frequency [Hz]
    configuration.minFreq        = min([minFreq maxFreq]);  % Minimum output spectrum frequency
    configuration.maxFreq        = max([minFreq maxFreq]);	% Maximum output spectrum frequency
    configuration.numPoints      = numPoints;               % Length desired of the output spectrum 

    % Calculate radar variables
    configuration.radarWaveLength = c./radarFrequency;                                                  % Radar Wavelength (lambda) [m]
    configuration.K               = 4*pi./configuration.radarWaveLength;                                 % Radar Wavenumber (K)
    configuration.frequency       = linspace(configuration.minFreq,configuration.maxFreq,numPoints);
    configuration.omega           = 2*pi.*configuration.frequency;                                       % Radian frequency (rad)
    
    % Generate struct for different species
    configuration.species = {};
    
    configuration.species.Name_s = {};              % Name_s - The species name
    configuration.species.N_s    = [];              % N_s - The density of species [m^-3]
    configuration.species.T_s    = [];              % T_s - Tempreture of species [ºK]
    configuration.species.V_s    = [];              % V_s - The Doppler velocity [m/s]
    configuration.species.q_s    = [];              % q_s - The charge of species [Coulomb]
    configuration.species.m_s    = [];              % m_s - Mass of specie [AMU]
    
    configuration.species.debyeLength_s  = [];      % Debye Length of Specie (lambda_D_s)
    configuration.species.thermalSpeed_s = [];      % Thermal Speed of Specie (C_s)
    configuration.species.omega_s        = [];      % Frequency of Specie (Omega_s)
    
    %Generate struct for electrons
    configuration.species.Name_e = 'E-';            % Name_e - The electron name
    configuration.species.N_e    = [];              % N_e - The density of species [m^-3]
    configuration.species.T_e    = [];              % T_e - Tempreture of species [ºK]
    configuration.species.V_e    = [];              % V_e - The Doppler velocity [m/s]
    configuration.species.q_e    = q_e;             % q_e - The charge of species [Coulomb]
    configuration.species.m_e    = m_e;             % m_E - Mass of specie [AMU]
    
    configuration.species.debyeLength_e  = [];      % Debye Length of electrons (lambda_De)
    configuration.species.thermalSpeed_e = [];      % Thermal Speed of electrons (C_e)
    configuration.species.omega_e        = [];      % Frequency of electrons (Omega_e)
   
    % Special options
    configuration.spectrum_output = 0;              % Spectrum_output - 0=all, 1=ion contribution, 2=electron contribution
end