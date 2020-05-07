
load short_modem_rx.mat

% The received signal includes a bunch of samples from before the
% transmission started so we need discard the samples from before
% the transmission started. 

start_idx = find_start_of_signal(y_r,x_sync);
% start_idx now contains the location in y_r where x_sync begins
% we need to offset by the length of x_sync to only include the signal
% we are interested in
y_t = y_r(start_idx+length(x_sync):end); % y_t is the signal which starts at the beginning of the transmission


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Put your decoder code here

%Plot of y_t 
plot(y_t)
%% Short Modem
t = [0:length(y_t)-1]'; 
y_c = y_t.*cos(2*pi*f_c*t/Fs);

filter_y = lowpass(y_c, .2);
subplot(2,2,1)
plot(y_t)
xlabel('Time (seconds)');
ylabel('Y_t')
subplot(2,2,2)
plot(t, y_c)
xlabel('Time (seconds)');
ylabel('Y_c');
subplot(2,2,[3,4])
plot(t, filter_y)
xlabel('Time (seconds)');
ylabel('Filtered Y_c');

%% 

data = filter_y %(50:4000,:)
T = 1/Fs;             % Sampling period       
L = length(data);     % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% n = length(data); 
% Ts = t(2)-t(1); % sample time
% Fs = 1/Ts; % sampling frequency
% NFFT = 2^nextpow2(n); % Next power of 2 from length of data
% Y = fft(data,NFFT)/n;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% Plot single-sided amplitude spectrum.
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

%%
vec_y = filter_y(50:100:4000,:);
vec_y_bits = (filter_y(50:100:4000,:)>=0);
t_vec = [0:length(vec_y)-1];
subplot(2,1,1)
stem(t_vec, vec_y)
xlabel('Time (seconds)');
ylabel('Filtered Y_c')
subplot(2,1,2)
stem(t_vec, vec_y_bits)
xlabel('Time (seconds)');
ylabel('Binary Filtered Y_c')


word = BitsToString(double(vec_y_bits))
%% Long modem
load long_modem_rx.mat

% The received signal includes a bunch of samples from before the
% transmission started so we need discard the samples from before
% the transmission started. 

start_idx = find_start_of_signal(y_r,x_sync);
% start_idx now contains the location in y_r where x_sync begins
% we need to offset by the length of x_sync to only include the signal
% we are interested in
y_t = y_r(start_idx+length(x_sync):end); % y_t is the signal which starts at the beginning of the transmission

%%
t = [0:length(y_t)-1]';
y_c = y_t.*cos(2*pi*f_c*t/Fs);
filter_y = lowpass(y_c, .2);
subplot(2,2,1)
plot(y_t)
xlabel('Time (seconds)');
ylabel('Y_t')
subplot(2,2,2)
plot(t, y_c)
xlabel('Time (seconds)');
ylabel('Y_c');
subplot(2,2,[3,4])
plot(t, filter_y)
xlabel('Time (seconds)');
ylabel('Filtered Y_c');
%%

vec_y = filter_y(50:100:(123*8*100),:)
vec_y_bits = (filter_y(50:100:(123*8*100),:)>=0)
t_vec = [0:length(vec_y)-1];
mystery = BitsToString(double(vec_y_bits))

subplot(2,1,1)
stem(t_vec, vec_y)
xlabel('Time (seconds)');
ylabel('Filtered Y_c')
subplot(2,1,2)
stem(t_vec, vec_y_bits)
xlabel('Time (seconds)');
ylabel('Binary Filtered Y_c')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% convert to a string assuming that x_d is a vector of 1s and 0s
% representing the decoded bits
BitsToString(x_d)

