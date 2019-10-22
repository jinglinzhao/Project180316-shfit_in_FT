% Use simulated spectra with planets AND stellar jitter: 
% /Volumes/DataSSD/SOAP_2/outputs/02.01/


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
K           = 2;
SN          = 10000;  
N_FILE      = 100;
t           = 1:N_FILE;
grid_size   = 0.1;
Fs          = 1/grid_size;
v0          = (-20 : grid_size : 20)';          % km/s

dir1        = '/Volumes/DataSSD/SOAP_2/outputs/02.01/';
dir2        = '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat/';
jitter      = importdata([dir1, 'RV.dat']) / 1000;      % activity induced RV [km/s]
% jitter      = jitter';
jitter      = [jitter', jitter'];
% jitter      = [jitter', jitter', jitter', jitter'];               % comment this out if not tesitng "planet + jitter"
idx         = (v0 >= -10) & (v0 <= 10);
v1          = v0(idx);


%%%%%%%%%%%%%%%%%%%
% Calculate Power %
%%%%%%%%%%%%%%%%%%%

% estimate the size of array FFT_power
filename    = [dir2, 'CCF', num2str(1), '.dat'];
A           = 1 - importdata(filename);
A           = A(idx);
A1          = A;
[~, bb, yy] = FUNCTION_FFT(A, Fs);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);
RV_noise    = zeros(1,N_FILE);
v_planet_array  = K * sin(t/100*4.5*2*pi + 1) * 0.001;     % comment this out if not tesitng "planet + jitter"
RV_gauss        = zeros(N_FILE,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked cross correlation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    v_planet    = v_planet_array(n) * 0;
    filename    = [dir2, 'CCF', num2str(mod(n-1,100)), '.dat'];
%     filename    = [dir2, 'CCF', num2str(1), '.dat'];        % choose the same line profile and shift it 
    A           = 1 - importdata(filename);
    A_spline    = spline(v0, A, v1-v_planet);
    % Plot the one without noise
    if mod(n,10) == 1
%         plot(v1, A_spline - A1, 'k')
        plot(v1, A_spline, 'k')
    end    
    % add noise
%     A_spline    = A_spline + normrnd(0, (1-A_spline).^0.5/SN);  
    
    % obtain line centroid 
    idx_fit     = (v1 >= -9) & (v1 <= 9);
    v_fit       = v1(idx_fit);
    A_spline_fit= A_spline(idx_fit);
    f           = fit(v_fit, A_spline_fit, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    RV_gauss(n) = f.b;
    
%     A_spline    = A_spline .* window;
    % NEW %
    if K == 10
        A_spline    = spline(v1, A_spline, v1+f.b);
    end
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);
    
end     
hold off
% title('Stacked cross correlation function')
% ylim([-1.2/1000 1.2/1000])
set(gca,'fontsize',20)
xlabel('Velocity [km/s]')
ylabel('Normalized intensity')
saveas(gcf,'1-Line_Profile','png')
% saveas(gcf,'SLPD1-Differential_line_Profile','png')
% saveas(gcf,'1-Differential_line_Profile','png')
% saveas(gcf,'LPD1-Line_Profile','png')
% saveas(gcf,'LPD1-Differential_line_Profile','png')
close(h)

% Determine the midpoint the equally divides the power spectrum %
cutoff_power= max(max(FFT_power)) * 0.0001; % used for publication for demonstration purpose
% cutoff_power= max(max(FFT_power)) * 0.001;
% cutoff_power= max(max(FFT_power)) * 0.0005;
% cutoff_power= max(max(FFT_power)) * 0.000005;
f_max       = max(FFT_frequency(FFT_power(:,1) > cutoff_power));
n           = abs(FFT_frequency) <= f_max;
power_sum   = sum(FFT_power(n,1));

% half power %
if 1
    cum = 0;
    for i = 1:fix(sum(n)/2)
        cum = cum + FFT_power(size(FFT_power,1)/2+1+i,1);
        if cum > power_sum/4
            break
        end
    end
    f_HL = FFT_frequency(size(FFT_power,1)/2+1+i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FT power in all epochs %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    h = figure; 
    hold on
    for n = 1:N_FILE
    %     plot(FFT_frequency, FFT_power(:, n) - FFT_power(:, 51), '-')
    %     title('Differential FT power in all epochs (overplot)')
    %     xlabel('FT frequency (1 / velocity in wavelength)')
    %     ylabel('Differential power') 

        plot(FFT_frequency, FFT_power(:, n), 'k')
    %     title('Stacked power spectrum')
    end 
    hold off
    xlabel('\xi [s/km]')
    ylabel('Power')   
    xlim([-0.15 0.1501])
    set(gca,'fontsize',20)
    saveas(gcf,'2-FT_power','png')
    % saveas(gcf,'2-Differential_FT_power','png')
    close(h)
end

if 1
    h = figure; 
    hold on
    for n = 1:N_FILE
        plot(FFT_frequency, FFT_power(:, n), 'k')
    end 
    grey = [0.4,0.4,0.4];
    plot([-f_max, -f_max], [0, 0.6], '--', 'Color', grey)
    plot([f_max, f_max], [0, 0.6], '--', 'Color', grey)
    plot([-f_HL, -f_HL], [0, 0.6], '--', 'Color', grey)
    plot([f_HL, f_HL], [0, 0.6], '--', 'Color', grey)
    hold off
    text(0, 0.3, 'L', 'FontSize', 20, 'HorizontalAlignment','center')
    text((f_max+f_HL)/2, 0.3, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
    text(-(f_max+f_HL)/2, 0.3, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
    text((f_max+0.2)/2, 0.3, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
    text(-(f_max+0.2)/2, 0.3, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
    xlabel('\xi [s/km]')
    ylabel('Power')   
    xlim([-0.199 0.199])
    set(gca,'fontsize',20)
%     saveas(gcf,'LPD2-FT_power','png')
    saveas(gcf,'2-FT_power','png')
    close(h)
end


%%%%%%%%%%%%%%%
% Phase angle %
%%%%%%%%%%%%%%%
h = figure;
plot(FFT_frequency, unwrap(angle(Y(:, 51))), '.')
title('Phase angle (Rotation phase = 0.51)')
xlabel('FT frequency (1 / velocity in wavelength)')
ylabel('Phase angle [radian]')
xlim([-0.15 0.15])
saveas(gcf,'3-Phase_angle','png')
close(h)


%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
n       = abs(FFT_frequency) <= f_max;
n2  = 0 <= FFT_frequency & FFT_frequency <= 0.1;
slope   = zeros(1,N_FILE);
RV_FT   = zeros(1,N_FILE);
RV_FT_err  = zeros(1,N_FILE);
Shift_FT = zeros(sum(n2), N_FILE);                                        % add this line @0411
% wegihted_velocity = zeros(1,N_FILE);
h = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(n);
    yy  = angle(Y(n, i)) - angle(Y(n, 1));
    xx2  = FFT_frequency(n2);
    yy2  = angle(Y(n2, i)) - angle(Y(n2, 1));    
    if mod(i,10) == 1
        plot(xx, yy, 'k-')
    end   
    % Phase angle -> RV
    weight          = FFT_power(n,i)';
    [fitresult, gof]= createFit(xx, yy', weight);
    slope(i)        = fitresult.p1;
    RV_FT(i)        = -slope(i) / (2*pi);
    Shift_FT(:,i)   = -gradient(yy2) ./ gradient(xx2') / (2*pi);              % add this line @0411
%     Shift_FT(:,i)   = -filter(-smooth_diff(2),1,yy2) ./ gradient(xx2') / (2*pi);              % add this line @0411
    ci              = confint(fitresult,0.95);
    RV_FT_err(i)    = abs(diff(ci(:,1))*1000 / (4*pi));
    if 0  % method 2
        idx_no0 = (xx~=0);
        velocity = - (yy(idx_no0)' ./ xx(idx_no0)) / (2*pi);
        wegihted_velocity(i) = sum(velocity .* weight(idx_no0)) ./ sum(weight(idx_no0)) * 1000;
    end
end
hold on 
% ymax = 0.0159;
% ymax = 0.009;
ymax = 0.011;
plot([-f_max, -f_max], [-ymax, ymax], '--', 'Color', grey)
plot([f_max, f_max], [-ymax, ymax], '--', 'Color', grey)
plot([-f_HL, -f_HL], [-ymax, ymax], '--', 'Color', grey)
plot([f_HL, f_HL], [-ymax, ymax], '--', 'Color', grey)
hold off
text(0, 0.005, 'L', 'FontSize', 20, 'HorizontalAlignment','center')
text((f_max+f_HL)/2, 0.005, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
text(-(f_max+f_HL)/2, -0.005, 'H', 'FontSize', 20, 'HorizontalAlignment','center')
text((f_max+0.2)/2, 0, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
text(-(f_max+0.2)/2, 0, '\oslash', 'FontSize', 20, 'HorizontalAlignment','center')
hold off
set(gca,'fontsize',20)
ylim([-ymax ymax])
xlim([-0.199 0.199])
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
saveas(gcf,'4-Relative_phase_angle','png')
% saveas(gcf,'LPD4-Relative_phase_angle','png')
% saveas(gcf,'SLPD4-Relative_phase_angle','png')
close(h)


%%%%%%%%%%%
% Plot 3D %
%%%%%%%%%%%
[X_fft,Y_fft] = meshgrid(t, FFT_frequency(n2));
h = surf(X_fft,Y_fft, Shift_FT*1000);
% ylim([0 0.1])
% set(h, 'edgecolor', 'none');
% h.FaceColor = 'interp';
% shading interp 
colorbar;
xlabel('JD');
xticks([]);
ylabel('\xi [s/km]');
zlabel('RV [m/s]');
% title('Line shift');
title('Line deformation');
set(gca,'fontsize',18)









