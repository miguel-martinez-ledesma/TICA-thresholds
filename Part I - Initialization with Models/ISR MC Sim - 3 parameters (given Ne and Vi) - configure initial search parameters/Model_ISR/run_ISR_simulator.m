% run_ISR_simulator 
%
%       Incoherent Scatter Radar Spectrum Simulator for MATLAB
% -------------------------------------------------------------------------
% Function to run the Incoherent Scatter Radar simulator for
% multiple Ion species
% -------------------------------------------------------------------------
% Inputs:
% -radarFreq Frequency of Radar (Hz)
% -NumPoints Number of points of spectrum
% -maxFreq Maximum frequency of spectrum
% -Te Tempreture of electrons[ºK]
% -Ti Tempreture of ions species [ºK]
% -Ne The density of electrons [m^-3]
% -Ve (optional) The electron Doppler velocity [m/s]
% -Vi (optional) The ion Doppler velocity [m/s]
% -p  (optional) Molecular Ion Fraction (p=O+/Ne)
%
% Outputs:
% -spectrum  Amplitude of ISR spectrum
% -frequency  Frequeny [Hz]
%
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 02 February 2016
% -------------------------------------------------------------------------


function [spectrum, frequency] = run_ISR_simulator(radarFreq,NumPoints,maxFreq,Te,Ti,Ne,Vi,p)

    % Check input arguments
    if (nargin > 8) || (nargin < 6)
        error('Wrong Number of Inputs: Function requires minimum 6 inputs (radarFreq, NumPoints, maxFreq, Te, Ti, Ne)');
    end
    % Fill in unset argument values.
    switch nargin
        case 6
            Vi = 0;
            p  = 0;
        case 7
            p  = 0;
    end

    % Calculate Ion density from Molecular Ion Fraction (p = [0 .. 1])
    % Assuming Sum(Ni)=Ne and O+ and "NO+ and O2+" (30.5[amu] as indicated in [Lathuillere et al., 1983]) 

    Ni_atomicOxigen  = Ne.*(1-p);
    Ni_molecularIons = Ne.*(p);

    % Nule Electron doppler velocity
    Ve = 0;

    %% Configuration

    % Initialize
    configuration = init_ISR_simulator(radarFreq, -maxFreq, maxFreq, NumPoints);

    % Add Ion Species
    configuration = addSpecie_ISR_simulator(configuration,'O+');
    configuration = addSpecie_ISR_simulator(configuration,'molecularIons');

    % Electrons configuration
    configuration = changeSpecie_ISR_simulator(configuration, 'E-', Ne, Te, Ve);

    % Ions configuration
    configuration = changeSpecie_ISR_simulator(configuration, 'O+', Ni_atomicOxigen, Ti, Vi);
    configuration = changeSpecie_ISR_simulator(configuration, 'molecularIons', Ni_molecularIons, Ti, Vi);

    %% ISR simulation
    
    spectrum  = ISR_simulator(configuration);
    frequency = configuration.frequency';
    