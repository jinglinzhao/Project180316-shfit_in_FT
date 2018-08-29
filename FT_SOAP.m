% Use simulated spectra with planets AND stellar jitter: 
% /Volumes/DataSSD/SOAP_2/outputs/02.01/


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN          = 10000;
N_FILE      = 100;
% N_FILE      = 200;
t           = 1:N_FILE;
grid_size   = 0.1;
Fs          = 1/grid_size;
v0          = (-20 : grid_size : 20)';          % km/s
dir1        = '/Volumes/DataSSD/SOAP_2/outputs/02.01/';
dir2        = '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat/';
% dir1      = '/Volumes/DataSSD/SOAP_2/outputs/HERMIT_2spot/';
% dir2      = '/Volumes/DataSSD/SOAP_2/outputs/HERMIT_2spot/fits/CCF_dat/';
jitter      = importdata([dir1, 'RV.dat']) / 1000;      % activity induced RV [km/s]
jitter      = jitter';
% jitter      = [jitter', jitter'];               % comment this out if not tesitng "planet + jitter"
idx         = (v0 >= -10) & (v0 <= 10);
v1          = v0(idx);

% window function %
if 0
    window  = v1 * 0 + 1;
    bound   = 8;
    idx_w   = abs(v1) >= bound;
    window(idx_w)   = (cos((abs(v1(idx_w))-bound)/(10-bound)*pi) + 1) /2;

    h = figure;
    plot(v1, window)
    title('Window function')
    xlabel('Wavelength in RV [km/s]')
    ylabel('Window function')
    saveas(gcf,'0-Window_function','png')
    close(h)
end

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
v_planet_array  = linspace(0,10,N_FILE) / 1000.;
% v_planet_array  = 2 * sin(t/100.*0.7*2*pi + 1) * 0.001;     % comment this out if not tesitng "planet + jitter"
RV_gauss        = zeros(N_FILE,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked cross correlation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    v_planet    = v_planet_array(n) * 0;
    filename    = [dir2, 'CCF', num2str(mod(n,100)), '.dat'];
%     filename    = [dir2, 'CCF', num2str(1), '.dat'];        % choose the same line profile and shift it 
    A           = 1 - importdata(filename);
    A_spline    = spline(v0, A, v1-v_planet);
    % Plot the one without noise
    if mod(n,10) == 1
        plot(v1, A_spline - A1, 'k')
    end    
    % add noise
    A_spline    = A_spline + normrnd(0, (1-A_spline).^0.5/SN);  
    
    % obtain line centroid 
    idx_fit     = (v1 >= -9) & (v1 <= 9);
    v_fit       = v1(idx_fit);
    A_spline_fit= A_spline(idx_fit);
    f           = fit(v_fit, A_spline_fit, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    RV_gauss(n) = f.b;
    
%     A_spline    = A_spline .* window;

    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);
    
end     
hold off
% title('Stacked cross correlation function')
set(gca,'fontsize',20)
xlabel('km/s')
ylabel('Normalized intensity')
% saveas(gcf,'1-Line_Profile','png')
% saveas(gcf,'1-Differential_line_Profile','png')
% saveas(gcf,'LPD1-Line_Profile','png')
saveas(gcf,'LPD1-Differential_line_Profile','png')
close(h)

% Determine the midpoint the equally divides the power spectrum %
f_max       = 0.15;
n           = abs(FFT_frequency) <= 0.15;
power_sum   = sum(FFT_power(n,1));
cum = 0;
for i = 0:fix(sum(n)/2)
    cum = cum + FFT_power(size(FFT_power,1)/2+i,1);
    if cum > power_sum/2
        break
    end
end
f_divide = FFT_frequency(size(FFT_power,1)/2+i);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FT power in all epochs %
%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% phase 2D plot 
if 0
    [aa,bb,cc] = FUNCTION_FFT(A1, Fs);
    h = figure;
    plot(real(cc), imag(cc), '.', real(Y(:,1)), imag(Y(:,1)), '.', 'MarkerSize', 10)
%     get(gca)
    legend('Noise free', 'SN = 10000')
    grid on 
    xlim([-0.002 0.008])
    ylim([-0.01 0.01])
    axis on
    xlabel('Real')    
    ylabel('Imaginary')
    title('Fourier transform in complex plane (zoomed-in)')
    set(gca,'fontsize',20)
    saveas(gcf,'7-Phase_angle_in_complex_plane_2','png')
%     close(h)
end


% Not used %
if 0
    figure; 
    hold on
    for n = (1024-25):(1024+25)
        i = 1:10;
        % plot(i, (angle(Y(n, :)) - angle(Y(n, 1))) / (FFT_frequency(n)+5.0049)  )
        plot(i, (angle(Y(n, i)) - angle(Y(n, 51))), '-')
        % plot(i, (angle(Y(n, :)) - angle(Y(n, 1))) / (  angle(Y(50, :)) - angle(Y(50, 1))  ), '.' )
        % disp( max(angle(Y(n, :))) - min(angle(Y(n, :))) )

        % plot(n, max(angle(Y(n, :))) - min(angle(Y(n, :))) , 'o')
    end 
    hold off
end


%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
n       = abs(FFT_frequency) <= f_max;
slope   = zeros(1,N_FILE);
RV_FT   = zeros(1,N_FILE);
RV_FT_err  = zeros(1,N_FILE);
% wegihted_velocity = zeros(1,N_FILE);
h = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(n);
    yy  = angle(Y(n, i)) - angle(Y(n, 1));
    if mod(i,10) == 1
        plot(xx, yy, 'k-')
    end
    % Phase angle -> RV
    weight          = FFT_power(n,i)';
    [fitresult, gof]= createFit(xx, yy', weight);
    slope(i)        = fitresult.p1;
    RV_FT(i)        = -slope(i) / (2*pi);
    ci              = confint(fitresult,0.95);
    RV_FT_err(i)    = abs(diff(ci(:,1))*1000 / (4*pi));
    if 0  % method 2
        idx_no0 = (xx~=0);
        velocity = - (yy(idx_no0)' ./ xx(idx_no0)) / (2*pi);
        wegihted_velocity(i) = sum(velocity .* weight(idx_no0)) ./ sum(weight(idx_no0)) * 1000;
    end
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
% saveas(gcf,'4-Relative_phase_angle','png')
saveas(gcf,'LPD4-Relative_phase_angle.png','png')
close(h)

% Low-pass %
nl      = (FFT_frequency >= 0) & (FFT_frequency < f_divide);
RV_FTL  = zeros(1,N_FILE);
RV_FTL_err  = zeros(1,N_FILE);
h       = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(nl);
    yy  = angle(Y(nl, i)) - angle(Y(nl, 1));
    if mod(i,10) == 1
        plot(xx, yy, 'k-')
    end
    % Phase angle -> RV
    weight          = FFT_power(nl,i)';

    if 0 % another way of fitting (from matlab manual but slower)
        ft = fittype( 'p1*x+p2', 'independent', 'x', 'dependent', 'y' );
        [xData, yData, weights] = prepareCurveData( xx, yy', weight);
        opts = fitoptions( ft );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf];
        opts.StartPoint = [0 0];
        opts.Upper = [Inf Inf];
        opts.Weights = weight;
        [fitresult, gof] = fit( xData, yData, ft, opts );
    end
    
    [fitresult, gof]= createFit(xx, yy', weight);
    slope(i)        = fitresult.p1;
    RV_FTL(i)       = -slope(i) / (2*pi);
    ci              = confint(fitresult,0.95);
    RV_FTL_err(i)   = abs(diff(ci(:,1))*1000 / (4*pi));    
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
title('Low-pass')
% saveas(gcf,'4-Relative_phase_angle_L','png')
saveas(gcf,'LPD4-Relative_phase_angle_L','png')
close(h)

% high-pass % 
n       = (FFT_frequency >= f_divide) & (FFT_frequency <= f_max);
RV_FTH  = zeros(1,N_FILE);
RV_FTH_err  = zeros(1,N_FILE);
h       = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(n);
    yy  = angle(Y(n, i)) - angle(Y(n, 1));
    if mod(i,10) == 1
        plot(xx, yy, 'k-')
    end
    % Phase angle -> RV
    weight          = FFT_power(n,i)';
    [fitresult, gof]= createFit(xx, yy', weight);
    slope(i)        = fitresult.p1;
    RV_FTH(i)       = -slope(i) / (2*pi);    
    ci              = confint(fitresult,0.95);
    RV_FTH_err(i)   = abs(diff(ci(:,1))*1000 / (4*pi));        
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
title('High-pass')
% saveas(gcf,'4-Relative_phase_angle_H','png')
saveas(gcf,'LPD4-Relative_phase_angle_H','png')
close(h)

% test
% figure; plot(1:100, wegihted_velocity, 1:100, RV_FT*1000)

%%%%%%%%%%%%%%%%%%%
% ONLY LINE SHIFT %
%%%%%%%%%%%%%%%%%%%
if 0
    % Compare with simulated RV % 
    h = figure;
        ax1 =subplot(20,1,1:17);
        xx = (v_planet_array-v_planet_array(1))*1000;
        yy1 = (RV_FT-RV_FT(1)) *1000;
        yy2 = (RV_gauss-RV_gauss(1))'*1000;
        [fitresult, gof]= createFit(xx, yy1, 1./(0.08+RV_FT_err).^2);
        fitresult
        hold on 
        scatter(xx, yy1, 'rs', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5)
        scatter(xx, yy2, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5)
        p0 = plot(xx, xx, 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
%         errorbar(xx, yy1, RV_FT_err, 'r.', 'MarkerSize', 0.1)
        hold off
        ylim([-0.1 10.1])
        title('RV Recovery')
        ylabel('Output RV [m/s]')
%         daspect(ax1,[1 1 1])
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'fontsize',15)
        legend({'RV_{FT}', 'RV_{Gaussian}', 'Output RV = Input RV'}, 'Location', 'northwest')

        positions = ax1.Position;
        ax2 = subplot(20,1,18:20);
        rms_gauss   = rms(yy2-xx - mean(yy2-xx));
        rms_FT      = rms(yy1-xx - mean(yy1-xx));  
        disp(rms_FT)
        hold on 
        scatter(xx, yy1-xx, 'rs', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5)
%         errorbar(xx, yy1-xx, RV_FT_err, 'r.', 'MarkerSize', 0.1)
        scatter(xx, yy2-xx, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5)
        p0 = plot([min(xx), max(xx)], [0,0], 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        xlabel('Input RV [m/s]')    
        ylabel('Residual [m/s]')
        set(gca,'fontsize',15)
        saveas(gcf,'5-LINE_SHIFT_ONLY','png')
    close(h)

    % high-pass and low-pass %
    h = figure;
        ax1 =subplot(20,1,1:17);
        xx = (v_planet_array-v_planet_array(1))*1000;
        yy1 = RV_FTH *1000;
        yy2 = RV_FTL *1000;
        [fitresult, gof]= createFit(xx, yy1, 1./(0.08+RV_FTH_err).^2);
        fitresult
        [fitresult, gof]= createFit(xx, yy2, 1./(0.08+RV_FTL_err).^2);
        fitresult        
        hold on 
        scatter(xx, yy1, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        scatter(xx, yy2, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        p0 = plot(xx, xx, 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        ylim([-0.1 10.1])
        title('RV Recovery')
        ylabel('Output RV [m/s]')
%         daspect(ax1,[1 1 1])
%         set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'fontsize',15)
        legend({'RV_{FT,H}', 'RV_{FT,L}', 'Output RV = Input RV'}, 'Location', 'northwest')

        positions = ax1.Position;
        ax2     = subplot(20,1,18:20);
        rms1    = rms(yy1-xx - mean(yy1-xx));  
        rms2    = rms(yy2-xx - mean(yy2-xx));
        hold on 
        scatter(xx, yy1-xx, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         errorbar(xx, yy1-xx, RV_FTH_err, 'k.', 'MarkerSize', 0.1)
        scatter(xx, yy2-xx, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         errorbar(xx, yy2-xx, RV_FTL_err, 'k.', 'MarkerSize', 0.1)
        p0 = plot([min(xx), max(xx)], [0,0], 'k--', 'LineWidth', 3); p0.Color(4)=0.3;
        hold off
        xlabel('Input RV [m/s]')    
        ylabel('Residual [m/s]')
        set(gca,'fontsize',15)
        saveas(gcf,'5-LINE_SHIFT_ONLY-HL','png')
    close(h)    
    
    
end    

%%%%%%%%%%%%%%%
% ONLY JITTER %
%%%%%%%%%%%%%%%
if 0
    % Compare with intrinsic line deformation RV % 
    h = figure; 
        yyL = RV_FTL*1000;
        yyH = RV_FTH*1000;
        xx = (RV_gauss - RV_gauss(1))'*1000;
        hold on
        p1 = scatter(xx, yyL, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
%         errorbar(xx, yyL, RV_FTL_err, 'k.', 'MarkerSize', 0.1)
        p1 = scatter(xx, yyH, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.4);
%         errorbar(xx, yyH, RV_FTH_err, 'k.', 'MarkerSize', 0.1)
        [fitresult_L, gof]= createFit(xx, yyL, 1./(1+RV_FTL_err*0).^2);
        fitresult_L        
        L1 = fitresult_L.p1;
        L2 = fitresult_L.p2;
        p3 = plot([min(xx), max(xx)], [L1*min(xx)+L2, L1*max(xx)+L2], 'k-', 'LineWidth', 2); p3.Color(4)=0.3;
        [fitresult_H, gof]= createFit(xx, yyH, 1./(1+RV_FTL_err).^2);
        fitresult_H        
        H1 = fitresult_H.p1;
        H2 = fitresult_H.p2;
        p4 = plot([min(xx), max(xx)], [H1*min(xx)+H2, H1*max(xx)+H2], 'k-', 'LineWidth', 2); p4.Color(4)=0.3;
        hold off
        xlabel('Jitter (RV_{Gaussian}) [m/s]')
        ylabel('RV_{FT} [m/s]')                     
        legend({'RV_{FT,L}', 'RV_{FT,H}'}, 'Location', 'northwest')
        set(gca,'fontsize', 15)
        saveas(gcf,'5-JITTER_ONLY_1','png')
    close(h)
    % p_fit =     0.7974   -0.0025 FOR A THREE-SPOT CONFIGRATION (0.8104)
    % p_fit =     0.7736   -0.2375 FOR A TWO-SPOT CONFIGRATION 

    if 0
        h = figure; 
            plot(xx, xx - yy, 'o')
            p_fit = polyfit(xx, xx - yy,1)
            title('FT with jitter only (residual)')
            xlabel('RV_{Gaussian} [m/s]')
            ylabel('RV_{Gaussian} - RV_{FT} (m/s)')
            saveas(gcf,'5-JITTER_ONLY_2','png')
        close(h)
    end
    
    % TIME SERIES
    h = figure;
    ax1 = subplot(3,1,1:2);
        xx = (1:N_FILE)/N_FILE;
        yy1 = RV_FT'*1000;
        yy2 = (RV_gauss - RV_gauss(1))*1000;
        hold on 
        scatter(xx, yy1, 'rs', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5);
        scatter(xx, yy2, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);
        hold off
        title('Apparent RV of deformed line profile')
        legend({'RV_{FT}', 'RV_{Gaussian}'}, 'Location', 'northwest')
        ylabel('RV [m/s]')
        set(gca,'fontsize',15)
    ax2 = subplot(3,1,3);
        scatter(xx, yy1 - yy2, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        xlabel('Stellar rotation phase')
        ylabel('RV difference [m/s]')
        set(gca,'fontsize',15)
    saveas(gcf,'5-JITTER_ONLY_3','png')
    close(h)
    
    % TIME SERIES 2 
    h = figure;
    ax1 = subplot(3,1,1:2);
        xx = (1:N_FILE)/N_FILE;
        yy = (RV_gauss - RV_gauss(1))*1000;
        hold on 
        scatter(xx, yy, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 1);
        scatter(xx, (yyL-L2)/L1, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(xx, (yyH-H2)/H1, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        hold off
        legend({'Jitter', 'RV_{FT,L} / k_{L}', 'RV_{FT,H} / k_{H}'}, 'Location', 'northwest')
        ylabel('RV [m/s]')
        set(gca,'fontsize',15)
    ax2 = subplot(3,1,3);
        hold on
        scatter(xx, (yyL-L2)/L1 - yy', 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(xx, (yyH-H2)/H1 - yy', 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        hold off
        rms((yyL-L2)/L1 - yy')
        rms((yyH-H2)/H1 - yy')
        xlabel('Stellar rotation phase')
        ylabel('Residual [m/s]')
        set(gca,'fontsize',15)
    saveas(gcf,'5-JITTER_ONLY_4','png')    
    close(h)    
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% LINE SHIFT AND JITTER %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare the total input RV and with the recovered RV
if 1

    % TIME SERIES
    h = figure; 
    ax1 = subplot(5,1,1:2);
        t    = (1:N_FILE)';
        yyL  = RV_FTL' * 1000;
        yyH  = RV_FTH' * 1000;
        yy2  = (RV_gauss - RV_gauss(1)) * 1000;
        hold on
        scatter(t, yyL, 15, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
%         scatter(t, yyH, 20, 'k*', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(t, yy2, 15, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);
        hold off
        title('RV recovery')
        ylabel('RV [m/s]')    
%         legend({'RV_{FT,L}', 'RV_{FT,H}', 'Gaussian'}, 'Location', 'north')
        legend({'RV_{FT,L}', 'Gaussian'}, 'Location', 'southwest')
        set(gca,'fontsize', 12)
        set(gca,'xticklabel',[])
%         dlmwrite('RV_IN.txt', yy2)
%         dlmwrite('RV_FT.txt', yy1)

        %     rv_g1 = sgolayfilt(rv_d,2,21);
        %     rv_g1 = sgolayfilt(rv_g1,2,11);
    ax2 = subplot(5,1,3:4); 
        rv_L        = yy2 - yyL;
        t_smooth    = linspace(1,200, 1000)';
        y_smooth1    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_L, t_smooth, 2);
        y_smooth11   = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_L, t, 2);
        xx2 = (jitter- jitter(1))' * 1000;
        p_fit1 = polyfit(xx2, rv_L, 1)
        hold on
        jitter_model1 = (rv_L-p_fit1(2))/p_fit1(1);
        c = @cmu.colors;
        plot(t, xx2, '--', 'color', [0.9100    0.4100    0.1700], 'LineWidth', 3)
        scatter(t, jitter_model1, 15, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)

%         rv_H        = yyH - yyL;
%         y_smooth2    = FUNCTION_GAUSSIAN_SMOOTHING(t, rv_H, t_smooth, 2);
%         p_fit2 = polyfit(xx2, rv_H, 1)
%         jitter_model2 = (rv_H-p_fit2(2))/p_fit2(1);
%         scatter(t, jitter_model2, 20, 'k+', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        
        plot1 = plot(t_smooth, (y_smooth1-p_fit1(2))/p_fit1(1), 'k', 'LineWidth', 2);
        plot1.Color(4) = 0.4;        
%         plot2 = plot(t_smooth, (y_smooth2-p_fit2(2))/p_fit2(1), 'k', 'LineWidth', 2);
%         plot2.Color(4) = 0.2;
%             pbaspect(ax2,[5 1 1])
        hold off
        ylabel('Jitter [m/s]')   
%         legend({'Jitter', 'Eq. 2.15', 'Eq. 2.16'}, 'Location', 'southwest')
        legend({'Jitter', 'Model'}, 'Location', 'northwest')
        set(gca,'fontsize', 12)
        set(gca,'xticklabel',[])

        ax3 = subplot(5,1,5);
        hold on 
        scatter(t, jitter_model1 - xx2, 15, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%         scatter(t, xx2 - jitter_model2, 20, 'k+', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
        plot1 = plot(t, (y_smooth11-p_fit1(2))/p_fit1(1) - xx2, 'k', 'LineWidth', 2);
        plot1.Color(4) = 0.4; 
        hold off
        xlabel('t')
        ylabel('Residual [m/s]')
        set(gca,'fontsize', 12)
        saveas(gcf,'5-PLANET_AND_JITTER2','png')

        rms(xx2 - mean(xx2))
        rms(jitter_model1 - xx2)
        rms((y_smooth11-p_fit1(2))/p_fit1(1) - xx2)
    close(h) 

    
    % Obtain kl, kh %
    yyL  = RV_FTL' * 1000;
    yyH  = RV_FTH' * 1000;
    yyG  = (RV_gauss - RV_gauss(1)) * 1000;
    plot(yyG-yyL, yyH-yyG, '.')
    p_fit = polyfit(yyG-yyL, yyH-yyG, 1)
    
    plot(yyG-yyL, yyH-yyL, '.')
    p_fit = polyfit(yyG-yyL, yyH-yyL, 1)    
    
    
    
    % REAL JITTER VS SCALED JITTER % 
    if 0
        h = figure; 
        hold on            
        plot(xx2, polyval(p_fit,xx2), '-')
        (sum((rv_d - polyval(p_fit,xx2)).^2)/N_FILE)^0.5
        (sum((rv_d - polyval(p_fit,xx2)).^2)/N_FILE)^0.5 / p_fit(1)
        hold off
        saveas(gcf,'7-Jitter_scaling','png')
        close(h) 
        % p_fit =    0.2019    0.0023 FOR A THREE-SPOT CONFIGRATION (0.1867)
        % p_fit =    0.2259    0.2377 FOR A TWO-SPOT CONFIGRATION
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    % NOT USED FOR NOW %
    figure; 
    xx = (jitter - v_planet_array(1:100)' - jitter(51))*1000;
    yy = slope/(-2*pi)*1000;
    plot(xx, yy' - xx, 'o')
    title('Recovered RV vs total RV')
    xlabel('input RV including jitter (m/s)')
    ylabel('residual (m/s)')    



    RVT_FT  = slope/(-2*pi)*1000;
    RVT     = (jitter - v_planet_array(1:100)' - jitter(51))*1000;

    figure; plot(t, RVT, t, RVT_FT')
    figure; plot(t, RVT - RVT_FT', t, (jitter-jitter(51))*1000)

    figure; plot(t, (RVT_FT' - RVT*0.7736)/(1-0.7736))




    x = (1:100)';
    y = (slope/(-2*pi)*1000)';
    rv = (jitter - v_planet_array(1:100)' - jitter(51))*1000;
    % dlmwrite('RV_tot.txt', rv)
    f = fit(x, y, 'a * sin(x/100*k+phi) * (1-m) + b + m * rv', 'StartPoint', [3 7 0.1 0.9 0 0]);






    ft = fittype('rv_ft(x, a, k, phi, m, b)');
    f = fit(t, RVT, ft);



    a0 = [3 7 0.1 0.9 0]; %?????a??????????
    options = statset('Jacobian', 'on');
    [a,r,J,~,msE] = nlinfit(x, y, @rv_ft, a0, options);%??
    [ypred,delta] = nlpredci(@rv_ft,x,a,r,'Jacobian',J,'predopt','observation','simopt','on');%??????????



    % test 

    x = (1:10)';
    y = x.^2+x/10;

    f = fit(x, y, 'a*x.^2+b', 'StartPoint', [0.01, 2]);
end



%%%%%%%%%%
% la fin %
%%%%%%%%%%

if 0


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct a template in FT space %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FFT_tpl = mean(FFT_power');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correct each observation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y_c = zeros(size(Y));
    Y_i = zeros(size(Y));
    iFFT_power = zeros(size(FFT_power));

    % Scale the power, and keep the phase

    for n = 1:N_FILE
        for i = 1:length(FFT_tpl)
            Y_c(i,n) = Y(i,n) / abs(Y(i,n)) *  FFT_tpl(i);
        end
    end

    for n = 1:N_FILE
        [iFFT_power(:, n), Y_i(:,n)] = FUNCTION_iFFT(Y_c(:,n));
    end

    figure;
    hold on
    for n = 1:N_FILE
        plot(FFT_frequency, iFFT_power(:, n) - iFFT_power(:, 1), '.')
    end
    hold off 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Examine the corrected line centre %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i_planet = zeros(1,N_FILE);
    for n = 1:N_FILE
        % idx     = iFFT_power(:, n) > 0.01;
        % f_power = fit( FFT_frequency(idx)'-0.5, iFFT_power(idx, n), 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [1, 0.5, 1, 0] );
        idx = (FFT_frequency<1) & (FFT_frequency>0);
        f_power = fit( FFT_frequency(idx)', iFFT_power(idx, n), 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0.5, 0.18, 0] );
        i_planet(n) = f_power.b;goo
    end

    % test %
    figure; plot(f_power, FFT_frequency(:), iFFT_power(:, n))

    v_planet    = 2 * sin(mod(MJD, 100)/100*7*2*pi + 1) * 0.001;
    figure; plot(RV_noise, i_planet, '.', 'markers',12)


end