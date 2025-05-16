
%% USER-DEFINED VARIABLES

w = -pi:(pi/100):pi;
% <-- Answer: Why is w from -pi to pi?
% 
% The variable 'w' represents the normalized radian frequency range for the frequency response analysis. 
% It is defined from -pi to pi to cover the full spectrum of frequencies for periodic signals, 
% allowing for the examination of both positive and negative frequencies, which is important in 
% signal processing to understand the behavior of filters across the entire frequency range.

%% HIGHPASS FILTER

% FREQUENCY RESPONSE
H2 = (1-exp(-1j*w*1));
% <-- Answer: What is the difference equation for this frequency response?

% -------------------> 
% The difference equation corresponding to this frequency response can be derived from the 
% transfer function H(z) = 1 - z^(-1). In the time domain, this can be expressed as
% y[n] = x[n] - x[n-1], where y[n] is the output and x[n] is the input signal. 
% This indicates that the output is the difference between the current input and the previous input, 
% allowing high-frequency components to pass while attenuating low-frequency components.
% <----------------------

% PLOT
figure;
subplot(2,1,1)
plot(w,abs(H2));    % ==> What does the abs() function do? 
                    % The abs() function computes the magnitude of the complex frequency response H2. 
                    % It returns the absolute value of each element in H2, which represents the amplitude of the 
                    % frequency response at each frequency 'w'. This is important for understanding how much 
                    % of each frequency component is passed through the filter.
                    % <==
grid on;
title('Magnitude Response')
xlabel('Normalized Radian Frequency');
ylabel('Amplitude');
subplot(2,1,2)
plot(w,angle(H2));  % ==> What does the angle() function do? 
                    % The angle() function calculates the phase angle (in radians) of the complex frequency response H2. 
                    % It returns the angle of each element in H2, which indicates the phase shift introduced by the 
                    % filter at each frequency 'w'. This is crucial for understanding how the filter affects the timing 
                    % of different frequency components in the input signal.
                    % <==
grid on;
title('Phase Response')
xlabel('Normalized Radian Frequency');
ylabel('Phase');

% <-- Answer: If you input a DC value into a highpass filter, what will be
%             its amplitude?
%
% The amplitude of a DC value (0 Hz frequency) input into a highpass filter will be 0. 
% This is because highpass filters are designed to attenuate low-frequency signals, including 
% DC components, effectively blocking them and allowing only higher frequency signals to pass through.