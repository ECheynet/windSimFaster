clearvars;close all;clc;
addpath ./functions
fs = 6; % sampling frequency (in Hertz)
M = 15; % M is the power of 2: If M = 10, 2^M = 1024 time step
[t,f] = getSamplingPara(M,fs);
N=numel(t); % number of time step (should be equal to 2^M)
Nfreq = numel(f); % number of frequency step (should be equal to 2^(M-1))
dt = median(diff(t)); % time step (should be equal to 1/fs)


Nyy = 5; % number of nodes along the y axis
Nzz = 5; % number of nodes along the z axis
ymin = -120;
ymax = 120;
zmin= 10;
zmax = zmin +100;
y = linspace(ymin,ymax,Nyy);
z = linspace(zmin,zmax,Nzz);
[Y,Z] = meshgrid(y,z);
M = numel(Z(:));


u_star = 0.8; % Friction velocity (m/s);
kappa = 0.4; % von karman constant
z0 = 0.01; % roughness length;
logProfile = @(u_star,z,z0,kappa)  u_star/kappa.*log(z./z0);
meanU = logProfile(u_star,Z,z0,kappa);

tic
Nf = numel(f);
u = zeros(N,M);
v = zeros(N,M);
w = zeros(N,M);
Su = zeros(Nf,M);
Sv = zeros(Nf,M);
Sw = zeros(Nf,M);
dummyU = meanU(:);
dummyZ = Z(:);

rng(1) % for reprodcuibility only!
for ii=1:M
    [Su(:,ii),Sv(:,ii),Sw(:,ii)] = KaimalModel(dummyU(ii),dummyZ(ii),f,u_star);
    [u(:,ii),t] = randomProcess(f,Su(:,ii));
    [v(:,ii),t] = randomProcess(f,Sv(:,ii));
    [w(:,ii),t] = randomProcess(f,Sw(:,ii));
end
toc

%%

Cuy = 10; % Davenport coherence coefficient for lateral separations and u component
Cuz = 12; % Davenport coherence coefficient for vertical separations and u component

Cvy = 10; % Davenport coherence coefficient for lateral separations and v component
Cvz = 12; % Davenport coherence coefficient for vertical separations and v component

Cwy = 5; % Davenport coherence coefficient for lateral separations and w component
Cwz = 6; % Davenport coherence coefficient for vertical separations and w component

[u_corr] = windSimFaster(Y,Z,meanU,Cuy,Cuz,f,u,'cohmodel','Davenport');
[v_corr] = windSimFaster(Y,Z,meanU,Cvy,Cvz,f,v,'cohmodel','Davenport');
[w_corr] = windSimFaster(Y,Z,meanU,Cwy,Cwz,f,w,'cohmodel','Davenport');
% reshape
u_corr = reshape(u_corr,[N,Nzz,Nyy]);
v_corr = reshape(v_corr,[N,Nzz,Nyy]);
w_corr = reshape(w_corr,[N,Nzz,Nyy]);
% Add the mean wind speed (only for the u-component)
u_corr = u_corr + reshape(meanU, [1, size(meanU)]);

if  Nyy==1,
    u_corr1(:,1,:)=u_corr;
    v_corr1(:,1,:)=v_corr;
    w_corr1(:,1,:)=w_corr;

    u_corr=u_corr1;
    v_corr=v_corr1;
    w_corr=w_corr1;   
end
%%


[~, indEnd] = min(abs(t-300));
clf;close all;
figure
tiledlayout(3,1,'TileSpacing','tight')

nexttile
plot(t(1:indEnd),squeeze(u_corr(1:indEnd,1,1:2:6)))
ylabel('u (m s^{-1})')
grid on; 
xlim([0 300])
ylim([-8 8]+min(meanU(:)))

nexttile
plot(t(1:indEnd),squeeze(v_corr(1:indEnd,1,1:2:6)))
ylabel('v (m s^{-1})')
grid on; 
xlim([0 300])
ylim([-8 8])

nexttile
plot(t(1:indEnd),squeeze(w_corr(1:indEnd,1,1:2:6)))
ylabel('w (m s^{-1})')
xlabel('time (s)')
grid on; 
xlim([0 300])
ylim([-8 8])

set(gcf, 'color', 'w'); % Set figure background color
set(findall(gcf, '-property',...
    'FontSize'), 'FontSize', 14, 'FontName', 'Times'); % Set font properties

%% Bowen

Cuy = [6 0];  % Davenport coherence coefficient for lateral separations and u component
Cuz = [6 18]; % Bowen coherence coefficient for vertical separations and u component

Cvy = [10 0]; % Davenport coherence coefficient for lateral separations and v component
Cvz = [6 18]; % Bowen coherence coefficient for vertical separations and v component

Cwy = [5 0];  % Davenport coherence coefficient for lateral separations and w component
Cwz = [3 18]; % Bowen coherence coefficient for vertical separations and w component

[u_corr] = windSimFaster(Y,Z,meanU,Cuy,Cuz,f,u,'cohmodel','Bowen');
[v_corr] = windSimFaster(Y,Z,meanU,Cvy,Cvz,f,v,'cohmodel','Bowen');
[w_corr] = windSimFaster(Y,Z,meanU,Cwy,Cwz,f,w,'cohmodel','Bowen');
% reshape
u_corr = reshape(u_corr,[N,Nzz,Nyy]);
v_corr = reshape(v_corr,[N,Nzz,Nyy]);
w_corr = reshape(w_corr,[N,Nzz,Nyy]);
% Add the mean wind speed (only for the u-component)
u_corr = u_corr + reshape(meanU, [1, size(meanU)]);

if  Nyy==1,
    u_corr1(:,1,:)=u_corr;
    v_corr1(:,1,:)=v_corr;
    w_corr1(:,1,:)=w_corr;

    u_corr=u_corr1;
    v_corr=v_corr1;
    w_corr=w_corr1;   
end

%%
[~, indEnd] = min(abs(t-300));
% clf;close all;
figure
tiledlayout(3,1,'TileSpacing','tight')

nexttile
plot(t(1:indEnd),squeeze(u_corr(1:indEnd,1,1:2:6)))
ylabel('u (m s^{-1})')
grid on; 
xlim([0 300])
ylim([-8 8]+min(meanU(:)))

nexttile
plot(t(1:indEnd),squeeze(v_corr(1:indEnd,1,1:2:6)))
ylabel('v (m s^{-1})')
grid on; 
xlim([0 300])
ylim([-8 8])

nexttile
plot(t(1:indEnd),squeeze(w_corr(1:indEnd,1,1:2:6)))
ylabel('w (m s^{-1})')
xlabel('time (s)')
grid on; 
xlim([0 300])
ylim([-8 8])

set(gcf, 'color', 'w'); % Set figure background color
set(findall(gcf, '-property',...
    'FontSize'), 'FontSize', 14, 'FontName', 'Times'); % Set font properties

%%

Nblock = 20; % a lot of block here to make it easier to visualize
Ncoh = round(N/Nblock);

distTarget = [25,50,75,100];
dy = abs(Y(:,1)'-Y(:,1)); % Matrix distance along y
dz = abs(Z(:,1)'-Z(:,1)); % Matrix distance along z
uCoh = 0.5*(meanU(:,1)'+meanU(:,1)); % Mean wind velocity between each nodes
z_avg = 0.5*(Z(:,1)'+Z(:,1));

%%%%%%%%%%%%%%%% Davenport %%%%%%%%%%%%%%%%%%%%%
% lateral separation
ay = Cuy(1).*dy;
% vertical separation
az = Cuz(1).*dz;
K_davenport = -sqrt(ay.^2+az.^2)./uCoh;
% Anonymous function

%%%%%%%%%%%%%%%%%% Bowen %%%%%%%%%%%%%%%%%%%%%%%%
% Bowen coefficient with lateral separation
ay = Cuy(1).*dy + Cuy(2).*dy.^2./z_avg;
% Bowen coefficient with lateral separation
az =Cuz(1).*dz + Cuz(2).*dz.^2./z_avg;
% Combine them into the coherence matrix for lateral and vertical separations
K_Bowen = -sqrt(ay.^2+az.^2)./uCoh;

cohFun = @(K,f) exp(K.*f);

Ntarget = numel(distTarget);

figure
tiledlayout(2,2,"TileSpacing","tight")
for ii=1:Ntarget
    [~,indDist] = min(abs(dz(1,:)-distTarget(ii)));
    [cocoh,~,freq] = coherence(u_corr(:,1,1),u_corr(:,indDist,1),Ncoh,round(Ncoh/2),Ncoh,fs);
    k = 2*pi*freq./uCoh(indDist,1);
    nexttile
    semilogx(k,cocoh,'ko','markerfacecolor','y');
    box on;hold on
    plot(k,cohFun(K_davenport(indDist,1),freq),'b','linewidth',1.5);
    plot(k,cohFun(K_Bowen(indDist,1),freq),'k','linewidth',1.5);
    legend(['dz = ',num2str(dz(indDist,1),2),' m'],'Davenport','Bowen','location','northeast');
    xlim([0,0.12])
    ylim([-0.3,1])
    ylabel('co-coherence')
    if ii>2,        xlabel('k (m^{-1})');    end
    set(gcf,'color','w')
end