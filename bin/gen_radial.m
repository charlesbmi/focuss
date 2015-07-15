function [k2, wi, varargout] = gen_radial(incr, N_r, N_proj, N_meas, shift)

%   If incr = 0, means Golden Ratio-based projection spacing
if incr == 0
    gr  =   0.5*(1+sqrt(5));
    incr=   180/gr;
end
if nargin == 4
    shift   =   0;
end

%   Generate angles
angles  =   0;
for i = 2:N_proj
    angles(i)   =   mod(angles(i-1)+incr,180);
end
angles = repmat(angles, [1, N_meas]);

%   Generate radial spacing information
r   =   linspace(-0.5,0.5,N_r);

%   Initialise output arrays
k   =   zeros(N_r, N_proj*N_meas);

%   Populate output arrays
for i = 1:N_proj*N_meas
   k(:,i)   =   (2*pi*r+shift)*exp(1j*pi*angles(i)/180);
end

%   Generate density compensation weighting 
%   Weight every point by |r| (dr, dthetha const) 
wi  =   abs(k);

k2  =   zeros(N_r, N_proj*N_meas, 2);
k2(:,:,1)   =   real(k);
k2(:,:,2)   =   imag(k);

if nargout == 3
    varargout{1}    =   angles;
end
