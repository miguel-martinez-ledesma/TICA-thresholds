% ISR_simulator
%
%       Incoherent Scatter Radar Spectrum Simulator for MATLAB
% -------------------------------------------------------------------------
% Info:
% Function to generate an Incoherent Scatter Radar spectrum for multiple
% Ion species implementing the Kudeki and Milla (2011) formulation for a
% non-collisional and non-magnetized plasma.
% -------------------------------------------------------------------------
% Based on:
%   -"Incoherent Scatter Radar Simulator (SimISR)" 
%       John Swoboda (Boston University)
%       https://github.com/jswoboda
% -------------------------------------------------------------------------
% References:
%   -"Incoherent Scatter Spectral Theories—Part I:
%       A General Framework and Results for Small Magnetic Aspect Angles" -
%       Erhan Kudeki and Marco A. Milla, 2011
% -------------------------------------------------------------------------
% -Outputs
%   Spectrum: Spectrum vector
%   freq    : Frequency vector [Hz]
%
% -Inputs
%   configuration : configurations
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 25 January 2016
% -------------------------------------------------------------------------

function [spectrum] = ISR_simulator(configuration)

    % Check number of different configurations
    num_config_s = size(configuration.species.N_s);
    num_config_e = length(configuration.species.N_e);
    num_config   = max([num_config_e num_config_s(2)]);

    for index = 1:num_config
        %--------------------------
        % ELECTRON

        % Compute Sigma_s and <|n_ts|^2>
        % Execute only once if one configuration is set
        if index <= length(configuration.species.N_e)
            % Get electron parameters
            N_s     = configuration.species.N_e(index);
            C_s     = configuration.species.thermalSpeed_e(index);
            omega_s = configuration.species.omega_e(index,:);
            debye_s = configuration.species.debyeLength_e(index);

            % Gordeyev Integral
            X_e       = sqrt(2) .* configuration.K .* C_s;
            theta_e   = omega_s ./ X_e;
            J_omega_e = ( sqrt(pi) .* exp(-theta_e .^ 2) - 1i.* 2 .* Faddeeva_Dawson(theta_e)) ./ X_e;

            % Conductance (Sigma_s)
            Sigma_e   = ((1j + omega_s .* J_omega_e) ./ ((configuration.K .^ 2) .* (debye_s .^ 2)));

            % Spectrum of the Thermal fluctuation of species <|n_ts|^2>
            n_te      = 2 .* N_s .* real(J_omega_e);
        end

        %--------------------------
        % IONS

        % Initialize
        omega_s     = zeros(1,configuration.numPoints);

        % Compute Sigma_s and <|n_ts|^2> for all ION species
        for s = 1:length(configuration.species.Name_s)
            % Execute only once if one configuration is set
            if index <= num_config_s(2)
                %Get ion parameters
                N_s        = configuration.species.N_s(s,index);
                C_s        = configuration.species.thermalSpeed_s(s,index);
                omega_s(:) = configuration.species.omega_s(s,index,:);
                debye_s    = configuration.species.debyeLength_s(s,index);

                % Gordeyev Integral
                X_s        = sqrt(2) .* configuration.K .* C_s;
                theta_s    = omega_s ./ X_s;
                J_omega_s  = ( sqrt(pi) .* exp(-theta_s .^ 2) - 1i.* 2 .* Faddeeva_Dawson(theta_s)) ./ X_s;

                % Conductance (Sigma_s)
                Sigma_s    = ((1j + omega_s .* J_omega_s) ./ ((configuration.K .^ 2) .* (debye_s .^ 2)));

                % Spectrum of the Thermal fluctuation of species <|n_ts|^2>
                n_ts       = 2 .* N_s .* real(J_omega_s);

                % Sum all sigmas and n_ts
                if (s==1)
                    Sum_Sigma_s = Sigma_s;
                    Sum_n_ts    = n_ts;
                else
                    Sum_Sigma_s = Sum_Sigma_s + Sigma_s;
                    Sum_n_ts    = Sum_n_ts + n_ts;
                end
            end
        end

        %--------------------------
        % Compute Spectrum
        denominator           = abs(1j + Sum_Sigma_s + Sigma_e) .^ 2;
        ion_contribution      = ((abs(Sigma_e) .^ 2) .* Sum_n_ts) ./ denominator;
        electron_contribution = ((abs(1j + Sum_Sigma_s) .^ 2) .* n_te) ./ denominator;

        if configuration.spectrum_output==0
            spectrum(:,index) = ion_contribution' + electron_contribution';
        elseif configuration.spectrum_output==1
            spectrum(:,index) = ion_contribution' ;
        elseif configuration.spectrum_output==2
            spectrum(:,index) = electron_contribution';
        end
        %--------------------------
    end
end