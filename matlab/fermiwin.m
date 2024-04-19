function w = fermiwin(N, r)
%% Fermi distribution window function
%  w = fermiwin(N, r)
%
% Returns a symmetric N-point Fermi window with transition region of width
% r*N/2.
%
% Input 
% -----
% N : positive integer
%   Window length (number of points)
% r : positive scalar
%   Transition region width factor.
%
% Output
% ------
% w : column vector
%   Fermi window
%
% See also TUKEYWIN, GAUSSWIN, WINDOW.
%
%% Created 2023-10-13 Samuel Adams-Tew

halfwidth = floor(N/2);
r_px = floor(halfwidth * r);

if r_px == 0
    w = ones(N, 1);
    return
end

% Create one half of the window
kT = r_px/9.2; % Speed of decay
ipt = halfwidth - 4.6*kT; % Inflection point
n = (1:halfwidth)';
w = 1./(1 + exp((n - ipt)/kT));

% Duplicate the window, and pad if needed
if halfwidth*2 == N
    w = [flipud(w); w];
else
    w = [flipud(w); 1; w];
end

end