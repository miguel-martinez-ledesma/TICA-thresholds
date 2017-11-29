% changeSpecie_ISR_simulator 
%
%       Incoherent Scatter Radar Spectrum Simulator for MATLAB
% -------------------------------------------------------------------------
% Function to change density, temperature, and doppler speed to species to
% the Incoherent Scatter Radar simulator for multiple Ion species
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 25 January 2016
% -------------------------------------------------------------------------



function [configuration] = changeSpecie_ISR_simulator(configuration, index, N_s, T_s, V_s)
    
    % Define generic constants
    kB   = 1.3807e-23;             % Boltzmann constant [J/ºK]
    eps0 = 8.854187817e-12;        % Permittivity of vacuum [F/m]

    % Verify if "index" is text
    if ~isnumeric(index)
        index_name = find(strcmpi(configuration.species.Name_s,index));
        if ~isempty(index_name)
            index = index_name;
        elseif strcmpi(index, 'E-')
            index=0;
        else
            error('ERROR: Unable to find Specie name: "%s" (Please add it first)\n',index)
        end
    end

    % Verify size of configurations
    num_config = max([length(N_s) length(T_s) length(V_s)]);
    if length(N_s) ~= num_config
        N_s = ones(1,num_config)*N_s;
    end
    if length(T_s) ~= num_config
        T_s = ones(1,num_config)*T_s;
    end
    if length(V_s) ~= num_config
        V_s = ones(1,num_config)*V_s;
    end

    if index == 0
        % ELECTRONS
        % Change electron values
        configuration.species.N_e = N_s;         % The density of species [m^-3]
        configuration.species.T_e = T_s;         % Tempreture of species [ºK]
        configuration.species.V_e = V_s;         % The Doppler velocity [m/s]

        % Calculate plasma characteristics
        % Debye Length of Specie (lambda_D_s)
        configuration.species.debyeLength_e = sqrt((eps0 .* kB .* configuration.species.T_e)./(N_s .* (configuration.species.q_e .^ 2)));
        % Thermal Speed of Specie (C_s)
        configuration.species.thermalSpeed_e = sqrt((kB .* T_s) ./ configuration.species.m_e);
        % Frequency of Specie (Omega_s)
        for i = 1:length(V_s)
            configuration.species.omega_e(i,:) = configuration.omega - configuration.K .* V_s(i);
        end
    else
        % IONS
        % Change specie values
        configuration.species.N_s(index,:) = N_s;         % The density of species [m^-3]
        configuration.species.T_s(index,:) = T_s;         % Tempreture of species [ºK]
        configuration.species.V_s(index,:) = V_s;         % The Doppler velocity [m/s]

        % Calculate plasma characteristics
        % Debye Length of Specie (lambda_D_s)
        configuration.species.debyeLength_s(index,:) = sqrt((eps0 .* kB .* T_s)./(N_s .* (configuration.species.q_s(index) .^ 2)));
        % Thermal Speed of Specie (C_s)
        configuration.species.thermalSpeed_s(index,:) = sqrt((kB .* T_s) ./ configuration.species.m_s(index));
        % Frequency of Specie (Omega_s)
        for i = 1:length(V_s)
            configuration.species.omega_s(index,i,:) = configuration.omega - configuration.K .* V_s(i);
        end
    end
end