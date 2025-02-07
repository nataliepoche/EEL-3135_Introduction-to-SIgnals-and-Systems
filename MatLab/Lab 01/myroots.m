function r = myroots(n, a)
% myroots: Find all the nth roots of the complex number a
%
% Input Args:
%   n: a positive integer specifying the nth roots
%   a: a complex number whose nth roots are to be returned
%
% Output:
%   r: 1xn vector containing all the nth roots of a 

% Converting to polar
A = abs(a);                                                     % Magnitude of coomplex number
phi = angle(a);                                                 % Argument (phase) of a

% nth root
roots = zeros(1,n);                                             % Initializes vector to stores the roots
for k = 0 : n-1;                                                % loops through n times

    % Formula for nth root
    r(k+1) = (A^(1/n)) * exp(1i * (phi + 2*pi*k)/n);            % r(k+1) is nth root (MATLAB starts from 1)
end

end