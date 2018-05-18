% Use simulated spectra with planets AND stellar jitter: 
% /Volumes/DataSSD/SOAP_2/outputs/02.01/


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN          = 10000;
N_FILE      = 200;
t           = 1:N_FILE;
grid_size   = 0.1;
Fs          = 1/grid_size;
v0          = (-20 : grid_size : 20)';          % km/s
dir1        = '/Volumes/DataSSD/SOAP_2/outputs/02.01/';
dir2        = '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat/';
% dir1      = '/Volumes/DataSSD/SOAP_2/outputs/HERMIT_2spot/';
% dir2      = '/Volumes/DataSSD/SOAP_2/outputs/HERMIT_2spot/fits/CCF_dat/';
jitter      = importdata([dir1, 'RV.dat']) / 1000;      % activity induced RV [km/s]
jitter      = [jitter', jitter'];
idx         = (v0 >= -10) & (v0 <= 10);
v1          = v0(idx);

% window function %
window  = v1 * 0 + 1;
bound   = 9;
idx_w   = abs(v1) >= bound;
window(idx_w)   = (cos((abs(v1(idx_w))-bound)/(10-bound)*pi) + 1) /2;

h = figure;
plot(v1, window)
title('Window function')
xlabel('Wavelength in RV [km/s]')
ylabel('Window function')
saveas(gcf,'0-Window_function','epsc')
close(h)
% window  = v1 * 0 + 1;

%%%%%%%%%%%%%%%%%%%
% Calculate Power %
%%%%%%%%%%%%%%%%%%%

% estimate the size of array FFT_power
filename    = [dir2, 'CCF', num2str(1), '.dat'];
A           = 1 - importdata(filename);
A           = A(idx);
A1          = A;
[aa, bb, yy]    = FUNCTION_FFT(A, Fs);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);
RV_noise    = zeros(1,N_FILE);
% v_planet_array  = linspace(-3,3,101) / 1000.;
v_planet_array  = 4 * sin(t/100.*1.8*2*pi + 1) * 0.001;
RV_gauss        = zeros(N_FILE,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked cross correlation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    v_planet    = v_planet_array(n);
    filename    = [dir2, 'CCF', num2str(mod(n,100)), '.dat'];
%     filename    = [dir2, 'CCF', num2str(1), '.dat'];        % choose the same line profile and shift it 
    A           = 1 - importdata(filename);
%     plot(v0(idx),A(idx)-A1, '.')
    A_spline    = spline(v0, A, v1-v_planet);
    
    % add noise
    A_spline    = A_spline + normrnd(0, (1-A_spline).^0.5/SN);  
    
    % obtain line centroid 
    idx_fit     = (v1 >= -9) & (v1 <= 9);
    v_fit       = v1(idx_fit);
    A_spline_fit= A_spline(idx_fit);
    f           = fit(v_fit, A_spline_fit, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    RV_gauss(n) = f.b;
    
    % FT
    A_spline    = A_spline .* window;
    plot(v1, A_spline)
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);
    
end     
hold off
title('Stacked cross correlation function')
xlabel('Wavelength in RV [km/s]')
ylabel('Normalized intensity')
saveas(gcf,'1-Line_Profile','epsc')
close(h)


%%%%%%%%%%%%%%%%%%%%%%%%%%
% FT power in all epochs %
%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure; 
hold on
for n = 1:N_FILE
    plot(FFT_frequency, FFT_power(:, n) - FFT_power(:, 51), '-')
    title('Differential FT power in all epochs (overplot)')
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('Differential power') 
    
%     plot(FFT_frequency, FFT_power(:, n))
%     title('FT power in all epochs (overplot)')
%     xlabel('FT frequency (1 / velocity in wavelength)')
%     ylabel('power')   
%     
    xlim([-0.2 0.2])
end 
hold off
saveas(gcf,'2-Differential_FT_power','epsc')
close(h)


%%%%%%%%%%%%%%%
% Phase angle %
%%%%%%%%%%%%%%%
h = figure;
plot(FFT_frequency, unwrap(angle(Y(:, 51))), '-')
title('Phase angle (Rotation phase = 0.51)')
xlabel('FT frequency (1 / velocity in wavelength)')
ylabel('Phase angle [radian]')
xlim([-0.35 0.35])
saveas(gcf,'3-Phase_angle','epsc')
close(h)


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
slope = zeros(1,N_FILE);
RV_FT = zeros(1,N_FILE);
h = figure; 
hold on
for i = 1:N_FILE
    % n = 1:50;
%     n = (1025-100):(1025+100);
%     n = (1025-40):(1025+40);    % plot for a particular frequency
%     n = (1025-10):(1025+10);          % THREE SPOT 
    n = (1025-5):(1025+5);          % TWO SPOT 
%     n = (size(FFT_frequency,2)/2-500):(size(FFT_frequency,2)/2+500);
    xx = FFT_frequency(n);
    yy = unwrap(angle(Y(n, i))) - unwrap(angle(Y(n, 1)));
    plot(xx, yy, '-')
    title('Phase angle (relative to 1st epoch)')
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('Phase angle [radian]')
    
    % Phase angle -> RV
    p = polyfit(xx, unwrap(angle(Y(n, i)))',1);
    slope(i) = p(1);
    RV_FT1 = -slope(1) / (2*pi);
    RV_FT(i) = -slope(i) / (2*pi) - RV_FT1;
end
hold off
saveas(gcf,'4-Relative_phase_angle','epsc')
close(h)


%%%%%%%%%%%%%%%%%%%
% ONLY LINE SHIFT %
%%%%%%%%%%%%%%%%%%%
if 0
    % Compare with simulated RV % 
    h = figure;
    xx = v_planet_array(1:100)*1000;
    xx = xx - xx(1);
    yy = RV_FT*1000;
    p_fit = polyfit(xx,yy,1)
    plot(xx, yy, 'o')
    title('Recovered RV vs input RV')
    xlabel('Input RV(m/s)')    
    ylabel('Recovered RV (m/s)')
    saveas(gcf,'5-LINE_SHIFT_ONLY_1','epsc')
    close(h)
    
    h = figure; 
    plot(xx, yy-xx, 'o')
    title('Recovered RV vs input RV (difference)')
    xlabel('Input RV (m/s)')    
    ylabel('Residual (m/s)')
    saveas(gcf,'5-LINE_SHIFT_ONLY_2','epsc')
    close(h)
end    
    
%%%%%%%%%%%%%%%
% ONLY JITTER %
%%%%%%%%%%%%%%%
if 0
    % Compare with intrinsic line deformation RV % 
    h = figure; 
    xx = (jitter - jitter(1))*1000;
    yy = RV_FT*1000;
    plot(xx, yy, 'o')
    title('FT with jitter only')
    xlabel('Jitter(m/s)')
    ylabel('RV_{FT} (m/s)')    
    p_fit = polyfit(xx,yy',1)
    saveas(gcf,'5-JITTER_ONLY_1','epsc')
    close(h)
    % p_fit =     0.7974   -0.0025 FOR A THREE-SPOT CONFIGRATION (0.8104)
    % p_fit =     0.7736   -0.2375 FOR A TWO-SPOT CONFIGRATION 

    h = figure; 
    plot(xx, yy-xx', 'o')
    title('FT with jitter only (residual)')
    xlabel('jitter(m/s)')
    ylabel('RV_{FT} - jitter (m/s)')
    saveas(gcf,'5-JITTER_ONLY_2','epsc')
    close(h)    
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% LINE SHIFT AND JITTER %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare the total input RV and with the recovered RV
h = figure; 
xx = (RV_gauss - RV_gauss(1)) * 1000;
yy = RV_FT*1000;
plot(xx, yy, 'o')
title('Planet and Jitter')
xlabel('Input RV including jitter (m/s)')
ylabel('RV_{FT} (m/s)')    
p_fit = polyfit(xx,yy',1)
dlmwrite('RV_IN.txt', xx)
dlmwrite('RV_FT.txt', yy')
saveas(gcf,'5-PLANET_AND_JITTER','epsc')
close(h)    


% TIME SERIES
h = figure; 
hold on
plot(t, xx'-yy, 'b')
xlabel('time')
ylabel('Scaled jitter (m/s)')   
plot(t, (jitter- jitter(1)) * 1000, 'r')
legend('(1-m)*jitter','jitter')
hold off
saveas(gcf,'6-Jitter_model','epsc')
close(h) 


% REAL JITTER VS SCALED JITTER % 
h = figure; 
xx2 = (jitter- jitter(1)) * 1000;
yy2 = xx'-yy;
plot(xx2, yy2, 'o')
xlabel('Real jitter (m/s)')
ylabel('Scaled jitter (m/s)')   
p_fit = polyfit(xx2,yy2, 1)
saveas(gcf,'7-Jitter_scaling','epsc')
close(h) 
% p_fit =    0.2019    0.0023 FOR A THREE-SPOT CONFIGRATION (0.1867)
% p_fit =    0.2259    0.2377 FOR A TWO-SPOT CONFIGRATION





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
        i_planet(n) = f_power.b;
    end

    % test %
    figure; plot(f_power, FFT_frequency(:), iFFT_power(:, n))

    v_planet    = 2 * sin(mod(MJD, 100)/100*7*2*pi + 1) * 0.001;
    figure; plot(RV_noise, i_planet, '.', 'markers',12)


end