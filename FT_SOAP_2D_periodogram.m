% 2D periodogram for the power

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 2000;
N_FILE          = 75;                               % number of CCF files
grid_size       = 0.1;
Fs          = 1/grid_size;
v0               = (-20 : grid_size : 20)';          % km/s
RV              = importdata('/Volumes/DataSSD/SOAP_2_copy/outputs/HERMIT_2spot_0720/RV.dat') / 1000;      % activity induced RV [km/s]
RV_gauss        = zeros(N_FILE,1);

idx             = (v0 > -10) & (v0 < 10);
v1               = v0(idx);

% estimate the size of array FFT_power
dir = '/Volumes/DataSSD/SOAP_2_copy/outputs/HERMIT_2spot_0720/CCF_dat/CCF';
filename    = [dir, num2str(1), '.dat'];
A           = 1 - importdata(filename);
A           = A(idx);
A1          = A;
[~, bb, yy] = FUNCTION_FFT(A, Fs);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);


% template %
A_tpl           = 1 - importdata('/Volumes/DataSSD/SOAP_2_copy/outputs/HERMIT_2spot_0720/CCF_tpl.dat');
A_tpl           = A_tpl(idx);
f_tpl           = fit( v, A_tpl, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
b_tpl           = f_tpl.b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked cross correlation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    i           = n - 1;
    v_planet    = 3 * sin(i/25*0.618*2*pi + 1) * 0.001 * 10; % km/s
    
    filename    = [dir, num2str(i), '.dat'];
%     filename    = [dir, num2str(1), '.dat'];        % choose the same line profile and shift it 
    A           = 1 - importdata(filename);
    A_spline    = spline(v0, A, v1-v_planet);
    % Plot the one without noise
    if mod(n,10) == 1
%         plot(v1, A_spline - A1, 'k')
        plot(v1, A_spline, 'k')
    end    
    % add noise
    A_spline    = A_spline + normrnd(0, (1-A_spline).^0.5/SN);  
    
    % obtain line centroid 
    idx_fit     = (v1 >= -9) & (v1 <= 9);
    v_fit       = v1(idx_fit);
    A_spline_fit= A_spline(idx_fit);
    f           = fit(v_fit, A_spline_fit, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    RV_gauss(n) = f.b;
    
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);
    
end     
hold off
% title('Stacked cross correlation function')
% ylim([-1.2/1000 1.2/1000])
set(gca,'fontsize',20)
xlabel('Velocity [km/s]')
ylabel('Normalized intensity')
% saveas(gcf,'1-Line_Profile','png')
close(h)

RV_gauss = RV_gauss - b_tpl;



array = 1:N_FILE;
[pxx_noise,f_noise] = plomb(FFT_power(0+1, array), array, 0.2, 100, 'normalized');
idx  = (FFT_frequency < 0.15) & (FFT_frequency >= 0);
% power = zeros(sum(idx), size(f_noise,1));
power = zeros(length(FFT_frequency), size(f_noise,1));

% simplified version % 
for i = 1:length(FFT_frequency)
    array = 1:N_FILE;
    [pxx_noise, f_noise] = plomb(FFT_power(i, array), array, 0.2, 100, 'normalized');
    power(i, :) = pxx_noise;
end

[X_fft,Y_fft] = meshgrid(f_noise, FFT_frequency(idx));
h = pcolor(X_fft,Y_fft, power(idx,:));
% ylim([0 0.1])
% zlim([min(min(Shift_FT)) max(max(Shift_FT))])
set(h, 'edgecolor', 'none');
h.FaceColor = 'interp';
colorbar;
xlabel('t');
ylabel('FT frequency domain');





