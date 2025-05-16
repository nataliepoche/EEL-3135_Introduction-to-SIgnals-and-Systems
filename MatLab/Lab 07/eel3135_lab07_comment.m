%% QUESTION #1 COMMENTING
clear
close all
clc

%% DEFINE FILTER AND INPUT
N = 100;
n = 0:(N-1);

% FILTER
b = (1/4)*[1 1 1]; % Numerator Coefficients of filter's transfer function, which is the coefficients of the input sequence in the difference equation
a = [1 -1 0.5 0.25]; % Denominator coefficients of filter's transfer function, which is the coefficients of output sequence in difference equation
% <-- Answer: What is the numerator of this filter's transfer function?
% The numerator of the filter's transfer function is (1/4) * [1 1 1] or
% (1/4) + (1/4)*z^-1 + (1/4)*z^-2

% <-- Answer: What is the denominator of this filter's transfer function?
% The denominator fo the filter transfer function is [1 -1 0.5 0.25] or 1 -
% z^-1 + 0.5*z^-2 + 0.25*z^-3

% Hint: look a later when we implement the filter with the "filter" command
%  (use "help filter" to see what it does)
%

% INPUT 1
x1 = zeros(N,1); 
x1(1) = 1; % Creates impulse input at n = 0
% <-- Answer: (True or False) x1 is an impulse input? If false,
% describe the input.
% True. x1 is an impulse input becuase it is a sequence of zeros with a
% single one at the first sample, definition of impulse

% INPUT 2
x2 = zeros(N,1); 
x2(1:12)  = cos(3*pi/2*   n(1:12)); % First segment of x2
x2(13:24) = cos(pi/4*   n(13:24)); % Second segment of x2
x2(25:36) = cos(2*pi*   n(25:36)); % Third segment of x2
x2(37:48) = cos(3*pi/4* n(37:48)); % Fourth segment of x2
x2(49:60) = cos(pi*     n(49:60)); % Fifth segment of x2
x2(61:72) = cos(pi/8* n(61:72)); % Sixth segment of x2
x2(73:84) = cos(pi/2*   n(73:84)); % Seventh segment of x2
% <-- Answer: (True or False) x2 is an single frequency input? If 
% false, describe the input.
% False. x2 is not a signle frequency input, it is a cobination of multiple
% cosine signals at different frequencies over time where each segment of
% x2 has a different frequency.

%% DEFINE AND PLOT OUTPUT

% OUTPUT 1
y1 = filter(b,a,x1);

% OUTPUT 2
y2 = filter(b,a,x2);


% PLOT THE IMPULSE RESPONSE AND DTFT
figure(1)
subplot(311)
stem(n,x1)
xlabel('Time (Samples)')
ylabel('x_1[n]')
subplot(312)
stem(n,y1)
xlabel('Time (Samples)')
ylabel('y_1[n]')
subplot(313)
pzplot(b,a)
axis equal

% PLOT THE IMPULSE RESPONSE AND DTFT
figure(2)
subplot(311)
stem(n,x2)
xlabel('Time (Samples)')
ylabel('x_2[n]')
subplot(312)
stem(n,y2)
xlabel('Time (Samples)')
ylabel('y_2[n]')
subplot(313)
pzplot(b,a)
axis equal



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