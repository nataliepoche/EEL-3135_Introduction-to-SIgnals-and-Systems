%% QUESTION #1 COMMENTING
clear
close all
clc

%% GET MUSIC

% The audio clip with from Xylo-Ziko's Album titled "Minimal" 
% from freemusicarchive.com. The Track name is "Phase 2"
% https://freemusicarchive.org/music/Xylo-Ziko/Minimal_1625/Phase_2_1586
[x, fs] = audioread('xylo-ziko.wav');


%% DEFINE TIME AND FREQUENCY AXES
w = -pi:pi/8000:pi-pi/8000;
t = 1/fs:1/fs:length(x)/fs;

%% DEFINE FILTER

% DEFINE POLES

mypoles = [ ...
        0.8*exp(1j*0.5e4/5e4*2*pi*0.9) ...
        0.8*exp(-1j*0.5e4/5e4*2*pi*0.9) ...        
        0.85*exp(1j*0.5e4/5e4*2*pi*0.95) ...
        0.85*exp(-1j*0.5e4/5e4*2*pi*0.95) ...
        0.85*exp(1j*0.5e4/5e4*2*pi*1.05) ...
        0.85*exp(-1j*0.5e4/5e4*2*pi*1.05) ...
        0.8*exp(1j*0.5e4/5e4*2*pi*1.1) ...
        0.8*exp(-1j*0.5e4/5e4*2*pi*1.1) ...
        0.9*exp(1j*0.5e4/5e4*2*pi) ...
        0.9*exp(-1j*0.5e4/5e4*2*pi) ...
         ... 
         ];
     
% DEFINE ZEROS     
myzeros = [ ...
        1 ...
        0.99*exp(1j*3726/5e4*2*pi) ...
        0.99*exp(-1j*3726/5e4*2*pi) ...
        0.99*exp(1j*5588/5e4*2*pi) ...
        0.99*exp(-1j*5588/5e4*2*pi) ...
        1*exp(1j*7471/5e4*2*pi) ...
        1*exp(-1j*7471/5e4*2*pi) ...
        0.96*exp(1j*9977/5e4*2*pi) ...
        0.96*exp(-1j*9977/5e4*2*pi) ...
        0.90*exp(1j*14960/5e4*2*pi) ...
        0.90*exp(-1j*14960/5e4*2*pi) ...
        0.85*exp(1j*9977/5e4*2*pi*2) ...
        0.85*exp(-1j*9977/5e4*2*pi*2) ...
         ... 
         ];


% CONVERT POLES AND ZEROS INTO COEFFICIENTS
[b,a] = pz2ba(mypoles,myzeros);

% COMPUTE GAIN TO MAINTAIN SIGNAL AMPLITUDE AROUND SOME FREQUENCY
G = abs(sum(a.*exp(-1j.*0.5e4/5e4*2*pi.*(0:(length(a)-1))))./sum(b.*exp(-1j.*0.5e4/5e4*2*pi.*(0:(length(b)-1)))));

% FILTER INPUT
y = filter(G*b,a,x);



%% PLOT FREQUENCY DOMAIN AND TIME DOMAIN INPUT / OUTPUT

% COMPUTE FREQUENY DOMAIN REPRESENTATION
X = DTFT(x,w);
Y = DTFT(y,w);

% PLOT FREQEUNCY DOMAIN AND
figure(1)
subplot(211)
plot(w,abs(X))
xlabel('Normalized Frequency (Radians)')
ylabel('Magnitude')
subplot(212)
plot(w,abs(Y))
xlabel('Normalized Frequency (Radians)')
ylabel('Magnitude')

% PLOT POLE-ZERO PLOT
figure(2)
pzplot(b,a)
axis equal;

% PLOT TIME DOMAIN
figure(3)
subplot(211)
plot(t,x)
xlabel('Time [s]')
ylabel('x[n]')
subplot(212)
plot(t,y)
xlabel('Time [s]')
ylabel('y[n]')

%% PLAY MUSIC

disp('Playing Original Music ... ')
soundsc(x, fs)
pause(length(x)/fs*1.1)

disp('Playing Filtered Music ... ')
soundsc(y, fs)

%% SAVE RESULT

% NORMALIZE THE SIGNAL TO AVOID CLIPPING
audiowrite('output.wav', y./max(abs(y)), fs)


%% ALL FUNCTIONS SUPPORTING THIS CODE %%
% ==================================================================
% NOTE: YOU DO NOT NEED TO ADD COMMENTS IN THE CODE BELOW. WE JUST 
% NEEDED POLE-ZERO PLOTTING CODE AND THUS WROTE IT. 
% ==================================================================
function pzplot(b,a)
% PZPLOT(B,A)  plots the pole-zero plot for the filter described by
% vectors A and B.  The filter is a "Direct Form II Transposed"
% implementation of the standard difference equation:
% 
%    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 

    % MODIFY THE POLYNOMIALS TO FIND THE ROOTS 
    b  = b(1:find(b,1,'last'));
    a  = a(1:find(a,1,'last'));
    b1 = zeros(max(length(a),length(b)),1); % Need to add zeros to get the right roots
    a1 = zeros(max(length(a),length(b)),1); % Need to add zeros to get the right roots
    b1(1:length(b)) = b;    % New a with all values
    a1(1:length(a)) = a;    % New a with all values

    % FIND THE ROOTS OF EACH POLYNOMIAL AND PLOT THE LOCATIONS OF THE ROOTS
    h1 = plot(real(roots(a1)), imag(roots(a1)));
    hold on;
    h2 = plot(real(roots(b1)), imag(roots(b1)));
    hold off;

    % DRAW THE UNIT CIRCLE
    circle(0,0,1)
    
    % MAKE THE POLES AND ZEROS X's AND O's
    set(h1, 'LineStyle', 'none', 'Marker', 'x', 'MarkerFaceColor','none', 'linewidth', 1.5, 'markersize', 8); 
    set(h2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor','none', 'linewidth', 1.5, 'markersize', 8); 
    axis equal;
    
    % DRAW VERTICAL AND HORIZONTAL LINES
    xminmax = xlim();
    yminmax = ylim();
    line([xminmax(1) xminmax(2)],[0 0], 'linestyle', ':', 'linewidth', 0.5, 'color', [1 1 1]*.1)
    line([0 0],[yminmax(1) yminmax(2)], 'linestyle', ':', 'linewidth', 0.5, 'color', [1 1 1]*.1)
    
    % ADD LABELS AND TITLE
    xlabel('Real Part')
    ylabel('Imaginary Part')
    title('Pole-Zero Plot')
    
end

function circle(x,y,r)
% CIRCLE(X,Y,R)  draws a circle with horizontal center X, vertical center
% Y, and radius R. 
%
    
    % ANGLES TO DRAW
    ang=0:0.01:2*pi; 
    
    % DEFINE LOCATIONS OF CIRCLE
    xp=r*cos(ang);
    yp=r*sin(ang);
    
    % PLOT CIRCLE
    hold on;
    plot(x+xp,y+yp, ':', 'linewidth', 0.5, 'color', [1 1 1]*.1);
    hold off;
    
end

function H = DTFT(x,w)
% DTFT(X,W)  compute the Discrete-time Fourier Transform of signal X
% acroess frequencies defined by W. 

    H = zeros(length(w),1);
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.'*(nn-1));
    end
    
end

function [b,a] = pz2ba(p,z)
% PZ2BA(P,Z) 	Converts poles P and zeros Z to filter coefficients
%               B and A

    % CONVERT ROOTS (POLES AND ZEROS) INTO POLYNOMIALS
    b = poly(z);
    a = poly(p);

end
