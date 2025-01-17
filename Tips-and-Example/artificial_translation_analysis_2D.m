% function artificial_translation_analysis_2D
% 
% Plot displacements after running DIC on translated images
% 
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2016

clear;
close all;
clc;

% Pixel sizes. Set to 1 to plot results in pixels
pix_size = 1; % um

% Expected translations
u_expected = [0 2 4 8];
v_expected = [0 2 4 8];

% File name with DIC results
dataname = 'DIC_results.mat';
% Name to save data. Set to empty array [] to prevent saving
savename = 'computational_translation_results_d0=06.png';

%% --- COMPUTE MEAN STRAINS ---

% Load data
load(dataname);

% Number of correlations
n_correlations = size(u,3);

% Preallocate
E11_mean = zeros(1,n_correlations);
E22_mean = E11_mean;
E12_mean = E11_mean;

% Run loop
for k=1:n_correlations
    u_k = u(:,:,k);
    v_k = v(:,:,k);
    [M, N] = size(u_k);
    % Only consider data that had a full subset to correlate
    idxx=ceil(w0/d0):floor(M-w0/d0);
    idxy=ceil(w0/d0):floor(N-w0/d0);
    ua = u_k(idxx,idxy);
    va = v_k(idxx,idxy);
    % Get rid of displacements that are too large
    ua(abs(ua)>w0/2)=nan;
    va(abs(va)>w0/2)=nan;
    
    % Compute normal strains
    [dudx, dudy] = grad2d(ua*pix_size/(d0*pix_size));
    [dvdx, dvdy] = grad2d(va*pix_size/(d0*pix_size));
    E11_mean(k) = nanmean(dudx(:));
    E22_mean(k) = nanmean(dvdy(:));
    E12_mean(k) = nanmean(1/2*(dudy(:)+dvdx(:)));
end

%% --- COMPUTE MEAN DISPLACEMENTS ---

for k=1:n_correlations
    u_k = u(:,:,k);
    v_k = v(:,:,k);
    [M, N] = size(u_k);
    % Only consider data that had a full subset to correlate
    idxx=ceil(w0/d0):floor(M-w0/d0);
    idxy=ceil(w0/d0):floor(N-w0/d0);
    ua = u_k(idxx,idxy);
    va = v_k(idxx,idxy);
    % Get rid of displacements that are too large
    ua(abs(ua)>w0/2)=nan;
    va(abs(va)>w0/2)=nan;
    
    % Stretch displacements by computed value of strain
    [M, N]=size(ua);
    %     % x strains corresponds to cols (n)
    %     for n=1:N, ua(:,n,:) = ua(:,n,:) - n*d0*E11_mean(k); end
    %     % y strains correspond to rows (m)
    %     for m=1:M, va(m,:,:) = va(m,:,:) - m*d0*E22_mean(k); end
    
    % Mean values
    u_m(k) = mean(ua(~isnan(ua)));
    v_m(k) = mean(va(~isnan(va)));
    
    %     % Subtract off first mean value
    %     if i==1
    %         um1 = u_m(1);
    %         vm1 = v_m(1);
    %     end
    %     u_m(i) = u_m(i)-um1;
    %     v_m(i) = v_m(i)-vm1;
    
    % Standard deviations
    u_std(k) = std(ua(~isnan(ua)));
    v_std(k) = std(va(~isnan(va)));
    
end

h1 = make_fig([0.6 0.6 1.2 1.2]);

subplot(2,2,1)
hold on
% X Translations
u_measured = u_m(1:4);
u_measured_std = u_std(1:4);
hE = errorbar(u_expected*pix_size,u_measured*pix_size,u_measured_std*pix_size,'o-r','LineWidth',1);
% Y Translations
v_measured = [v_m(1) v_m(5:7)];
v_measured_std = [v_std(1) v_std(5:7)];
hE = errorbar(v_expected*pix_size,v_measured*pix_size,v_measured_std*pix_size,'^-b','LineWidth',1);
% Plot Settings
h = legend('U_x','U_y','location','northwest');
if pix_size==1
    xlabel('Applied displacement (pix)');
    ylabel('Measured displacement (pix)');
else
    xlabel('Applied displacement (\mum)');
    ylabel('Measured displacement (\mum)');
end
set(h,'box','off');
set(gca,'box','off');
% set(gca,'xtick',0:1:10,'ytick',0:1:10);
% axis([0 6 0 10])
xlim([0 10])
hold off

subplot(2,2,2) % ERRORS
u_err = abs(u_measured-u_expected)./u_expected;
v_err = abs(v_measured-v_expected)./v_expected;
hold on
plot(u_expected*pix_size,u_err*100,'o-r');
plot(v_expected*pix_size,v_err*100,'^-b');
if pix_size==1
    xlabel('Applied displacement (pix)');
else
    xlabel('Applied displacement (\mum)');
end
ylabel('Error (%)');
h = legend('U_x','U_y','location','northeast');
set(h,'box','off');
set(gca,'box','off');
xlim([0 10]);

% --- PLOT STRAINS ---

subplot(2,2,3);
hold on;
plot(E11_mean*100,'o-');
plot(E22_mean*100,'^-');
plot(E12_mean*100,'s-');
xlabel('Image no.')
ylabel('Strain (%)');
h = legend('\epsilon_{11}','\epsilon_{22}','\epsilon_{12}','location','southwest');
set(h,'box','off');

% Save plots
if ~isempty(savename)
    print('-dpng','-r300',savename);
    % print('-depsc',savename)
end
