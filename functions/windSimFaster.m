function u_corr=windSimFaster(Y,Z,U,Cy,Cz,f,u,Uhub,z_hub,varargin)
%WINDSIMFASTER Introduce spatial coherence into uncorrelated wind histories.
%
% u_corr=windSimFaster(Y,Z,U,Cy,Cz,f,u,Uhub,z_hub,'cohmodel',model)
% applies the selected coherence model to the columns of u. The same
% function can be called independently for the u, v, and w components.
%
% Inputs:
%   Y,Z     - Grid-coordinate matrices [m].
%   U       - Mean wind speed at each grid node [m/s].
%   Cy,Cz   - Lateral and vertical model coefficients. These inputs are not
%             used by the IEC model.
%   f       - Positive frequency vector [Hz].
%   u       - Uncorrelated histories [time samples x grid nodes].
%   Uhub    - Hub-height mean wind speed [m/s].
%   z_hub   - Hub height [m].
%   model   - 'Davenport', 'Vogt', or 'IEC'.
%
% Output:
%   u_corr  - Spatially correlated histories with the same size as u.
%
% For the IEC option, the IEC exponential co-coherence is applied without
% component-specific parameters. Call this function with the same Y and Z
% grids for u, v, and w to apply the same IEC coherence to all components.
%
% Author: E. Cheynet - UiB. Revised 2026-09-04.

%% Parse and validate inputs
p=inputParser;
p.CaseSensitive=false;
p.addParameter('cohmodel','Davenport',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
p.parse(varargin{:});
cohmodel=validatestring(p.Results.cohmodel,{'Davenport','Vogt','IEC'});

if numel(Y)~=numel(Z) || numel(Y)~=numel(U)
    error('windSimFaster:GridSizeMismatch','Y, Z, and U must contain the same number of grid nodes.')
end
if size(u,2)~=numel(Y)
    error('windSimFaster:HistorySizeMismatch','The number of columns in u must equal the number of grid nodes.')
end
if Uhub<=0 || z_hub<=0
    error('windSimFaster:InvalidHubParameters','Uhub and z_hub must be positive.')
end

%% Prepare pairwise grid quantities
dy=abs(Y(:)'-Y(:));
dz=abs(Z(:)'-Z(:));
z_avg=0.5.*(Z(:)'+Z(:));
meanU=0.5.*(U(:)'+U(:));

%% Select coherence model
if strcmpi(cohmodel,'Davenport')
    ay=Cy(1).*dy;
    az=Cz(1).*dz;
    K=-sqrt(ay.^2+az.^2)./meanU;
    modelFunCoh=@(K,freq)exp(K.*freq);
elseif strcmpi(cohmodel,'Vogt')
    ay=dy.*Cy(1).*exp(Cy(2).*dy./z_avg);
    az=dz.*Cz(1).*exp(Cz(2).*dz./z_avg);
    M=size(dy,1);
    K=nan(M,M,2);
    K(:,:,1)=-sqrt(ay.^2+az.^2)./meanU;
    K(:,:,2)=-sqrt(Cy(3).^2.*dy.^2+Cz(3).^2.*dz.^2)./z_avg;
    modelFunCoh=@(K,freq)exp(K(:,:,1).*freq+K(:,:,2));
elseif strcmpi(cohmodel,'IEC')
    K=sqrt(dy.^2+dz.^2);
    Lambda1=0.7.*min(60,z_hub);
    Lc=8.1.*Lambda1;
    modelFunCoh=@(K,freq)exp(-12.*sqrt((freq.*K./Uhub).^2+(0.12.*K./Lc).^2));
end

%% Introduce spatial coherence
N=size(u,1);
if mod(N,2)~=0
    error('windSimFaster:OddSampleCount','The number of time samples must be even.')
end
df=median(diff(f));
k=0:N-1;
f2s=min(k,N-k).*df;
fftU0=fft(u);
fftU=fftU0;
for ii=2:N/2
    cohU=modelFunCoh(K,f2s(ii));
    try
        C=chol(cohU,'lower');
    catch
        [L,D]=ldl(cohU,'lower');
        C=L*sqrt(D);
    end

    phaseVector = exp(1i.*angle(fftU0(ii,:)'));
    dummy = C*phaseVector;
    fftU(ii,:) = abs(fftU0(ii,:)).*dummy';
    fftU(N-ii+2,:) = conj(fftU(ii,:));


end
fftU(1,:)=0;
fftU(N/2+1,:)=real(fftU0(N/2+1,:));
u_corr=ifft(fftU,'symmetric');
end
