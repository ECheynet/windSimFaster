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
% Author: E. Cheynet - uiB - Last modified 2024-04-16
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
else
    error('unknown coherence model')
end




%% Introduce coherence
% Define a frequency vector
f2s = [0 f(1:end-1)];
f2s = [f2s fliplr(f2s)];% 2-sided frequency array
tic
N = numel(f2s);
fftU = (fft(u));
for ii=2:N
    cohU = modelFunCoh(K,f2s(ii));
    try
        C = chol(cohU,'lower');
    catch exception
        [L,D]=ldl(cohU,'lower'); % a LDL decomposition is applied this time
        C = (L*sqrt(D));
    end
    fftU(ii,:)= C*exp(1i*angle(fftU(ii,:)'));
end
fftU = fftU.*abs(fft(u));
toc
% zero mean value
fftU(1,:)=0;
% Get back to the time domain
u_corr = real(ifft(fftU));
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