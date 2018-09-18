%% Plot CSDs from all samples
clear;  close all


%% Load data and calculate distribution
CSD = load('constrictions_Zhejiang');
CSD = CSD.constrictions(:,4);
CSD = sortrows(CSD);


PSD = load('particles_Zhejiang');
PSD = PSD.particledata(:,5);
PSD = sortrows(PSD);


CSD_dist = (100*(1:length(CSD))./length(CSD))';

vol = (PSD(:,1).^3) * (4/3) * (pi);
PSD_dist = 100 * cumsum(vol) / sum(vol);


%% Option to normalise by smallest particles
CSD  = CSD(:,1) ./ min(PSD);
PSD = PSD ./ min(PSD)


%%
figure;



hold on; box on;
xlabel('Radius (m)')
ylabel('% smaller')

plot(CSD, CSD_dist, 'color', 'r', 'linewidth', 3);
plot(PSD, PSD_dist, 'color', 'k', 'linewidth', 3);


legend('CSD', 'PSD')
set(gca, 'xscale', 'log')
% xlim([0.155 0.8])



return;

