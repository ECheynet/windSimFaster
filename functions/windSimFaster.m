function [u_corr] = windSimFaster(Y,Z,U,Cy,Cz,f,u,varargin)
% windSimFaster Introduces spatial correlation into uncorrelated wind
%velocity histories using Davenport's coherence model.
%
%   This function simulates correlated wind fields over a grid defined by Y
%   and Z coordinates, based on mean wind speeds U, Davenport coherence
%   coefficients Cy and Cz, frequency vector f, and an array of uncorrelated
%   wind velocities u.
%
% Syntax:
%   [u_corr] = windSimFaster(Y, Z, U, Cy, Cz, f, u)
%
% Inputs:
%   Y   - Matrix of y-coordinates (lateral) of grid nodes. [Nzz,Nyy] = size(Y).
%   Z   - Matrix of z-coordinates (vertical) of grid nodes. Size must match Y.
%   U   - Vector or matrix of mean wind speed at each grid node. If a vector,
%         it must have one element for each node.
%   Cy  - Scalar Davenport coherence coefficient for lateral (y-direction)
%         separations.
%   Cz  - Scalar Davenport coherence coefficient for vertical (z-direction)
%         separations.
%   f   - Vector of frequency steps for the simulation.
%   u   - Matrix of uncorrelated wind velocities. Each column corresponds to
%         a grid node and each row to a time step.
%
% Outputs:
%   u_corr - Matrix of correlated wind velocities. Same dimensions as input
%            'u', with each column representing a grid node and each row
%            representing a time step.
%
% Example:
% % Define grid coordinates
% Y = linspace(-50, 50, 20);
% Z = linspace(5, 100, 20);
% [Y, Z] = meshgrid(Y, Z);
% % Mean wind speed at each grid node
% U = 10*ones(size(Y));
% % Davenport coherence coefficients
% Cy = 3;
% Cz = 6;
% % Frequency vector
% f = linspace(0, 1, 100);
% % Uncorrelated wind velocities
% u = randn(2*length(f), numel(Y));
% % Simulate spatially correlated wind velocities
% u_corr = windSimFaster(Y, Z, U, Cy, Cz, f, u);
%
% Author: E. Cheynet - uiB - Last modified 2026-07-23
%
% See also FFT, IFFT, CHOL, LDL.

%% input parser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('cohmodel','Davenport');
p.parse(varargin{:});
cohmodel = p.Results.cohmodel;
%% Prepare data for the coherence
% [Nzz,Nyy]=size(Y);
dy = abs(Y(:)'-Y(:)); % Matrix distance along y
dz = abs(Z(:)'-Z(:)); % Matrix distance along z
z_avg = 0.5*(Z(:)'+Z(:));
meanU = 0.5*(U(:)'+U(:)); % Mean wind velocity between each nodes

% Anonymous function
if strcmpi(cohmodel,'Davenport')

    % Davenport decay coefficient Cy(1x1) with lateral separation
    ay = Cy(1).*dy;
    % Davenport decay coefficient Cz (1x1) with  vertical separation
    az = Cz(1).*dz;
    % Combine them into the coherence matrix for lateral and vertical
    % separations
    K = -sqrt(ay.^2+az.^2)./meanU;

    modelFunCoh = @(K,f) exp(K.*f);

elseif strcmpi(cohmodel,'Bowen')

    % Bowen coefficient with lateral separation
    ay = Cy(1).*dy + Cy(2).*dy.^2./z_avg;
    % Bowen coefficient with lateral separation
    az =Cz(1).*dz + Cz(2).*dz.^2./z_avg;

    % Combine them into the coherence matrix for lateral and vertical separations
    K = -sqrt(ay.^2+az.^2)./meanU;

    modelFunCoh = @(K,f) exp(K.*f);


elseif strcmpi(cohmodel,'Bowen2')
    [K] = coh_bowen2(Cy,Cz,f,meanU,dy,dz,z_avg);
    modelFunCoh = @(K,f) exp(K.*f);

elseif strcmpi(cohmodel,'vogt')

    % Mod. Bowen coefficient with lateral separation
    ay = dy.* Cy(1).* exp(Cy(2).*dy./z_avg);
    % Mod. Bowen coefficient with vertical separation
    az = dz.* Cz(1).* exp(Cz(2).*dz./z_avg);

    % Combine them into the coherence matrix for lateral and vertical separations
    M = numel(dy(:,1));
    K = nan(M,M,2);
    K(:,:,1) = -sqrt(ay.^2+az.^2)./meanU;

    dummy1 = Cy(3).^2*dy.^2;
    dummy2 = Cz(3).^2*dz.^2;
    K(:,:,2) = -1./z_avg.*sqrt(dummy1 + dummy2);

    modelFunCoh = @(K,f) exp(K(:,:,1).*f +K(:,:,2).*ones(size(f)));
elseif strcmpi(cohmodel,'IEC')
else
    error('unknown coherence model')
end




%% Introduce coherence
% Define a frequency vector
tic
N = size(u,1);
df = median(diff(f));
k = 0:N-1;
f2s = min(k,N-k)*df;
fftU0 = fft(u);
fftU = fftU0;
for ii = 2:N/2
    cohU = modelFunCoh(K,f2s(ii));
    C = chol(cohU,'lower');

    dummy = C*exp(1i*angle(fftU0(ii,:)'));
    fftU(ii,:) = abs(fftU0(ii,:)).*exp(1i*angle(dummy.'));
    fftU(N-ii+2,:) = conj(fftU(ii,:));
end
fftU(1,:) = 0;
fftU(N/2+1,:) = real(fftU0(N/2+1,:));
u_corr = ifft(fftU,'symmetric');


% % Add the mean wind speed
% u_corr = (u_corr'+ U(:))';



%% Nested function

    function [K_Bowen2] = coh_bowen2(Cy,Cz,f,meanU,dy,dz,Z)


        % Modified Bowen coefficient with lateral separation
        ay_Bowen2 = sqrt((Cy(1).*f).^2 + (Cy(3)).^2).*dy + Cy(2).*dy^2./Z.*f;
        % Bowen coefficient with lateral separation
        az_Bowen2 = sqrt((Cz(1).*f).^2 + (Cz(3)).^2).*dz + Cz(2).*dz^2./Z.*f;
        % Combine them into the coherence matrix for lateral and vertical separations
        K_Bowen2 = -sqrt(ay_Bowen2.^2+az_Bowen2.^2)./(f.*meanU);
    end
end
