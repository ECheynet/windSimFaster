function [Su,Sv,Sw] = KaimalModel(U,Z,f,u_star)
% [Su,Sv,Sw,Suw,Svw] = KaimalModel(U,Z,f,u_star) computes the one-point
% auto and cross-spectral densities of the Kaimal model [1].
% [1] Kaimal, J. C., & Finnigan, J. J. (1994).
% Atmospheric boundary layer flows: their structure and measurement.
% Oxford university press.
% 
% Inputs:
% U: scalar [1x1] of mean wind velocity (in m/s) at each node of a grid.
% Z: scalar [1x1] of height (in m) at each node of a grid.
% f: vector [1 x Nfreq] of frequency (in Hz)
% u_star: scalar [1 x 1] friction velocity (in m/s)
%
% Outputs:
% Recalling that PSD = power spectral density and
% CPSD = cross-power spectral density:
% Su: vector [1 x Nfreq] corresponding to the PSD the u-component
% Sv: vector [1 x Nfreq] corresponding to the PSD the v-component.
% Sw: vector [1 x Nfreq] corresponding to the PSD the w-component
% Suw: vector [1 x Nfreq] corresponding to the CPSD the u and w components
%
% Author: E. Cheynet - UiB - last modified : 10-10-2025

%%
fr = (f*Z/U);
scalingCoeff = u_star.^2./f;
Su = 102.*fr./(1+33.*fr).^(5/3).*scalingCoeff; % Kaimal  model (NOT normalized)
Sv = 17.*fr./(1+9.5.*fr).^(5/3).*scalingCoeff; % Kaimal  model (NOT normalized)
Sw = (2.*fr./(1+5.*fr.^(5/3))).*scalingCoeff; % Kaimal model (NOT normalized)
end