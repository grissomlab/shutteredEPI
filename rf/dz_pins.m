function [rf,g,ttipdown,Npulses,NtSub] = dz_pins(tb,slsep,dthick,de,maxb1,gmax,gslew,dt,minRFdur,halfShift,pType)

% design a PINS pulse
% 10/1/12 WA Grissom

% Output variables are: rf, g
%

if ~exist('halfShift','var')
  halfShift = false; % do not shift profile by half a slice gap
end

if ~exist('pType','var')
    pType = 'se';
end

kzw = tb/(dthick/10); % 1/cm, width in kz-space we must go
Npulses = ceil(kzw/(1/slsep)); % number of subpulses

% call SLR to get pulse
rfSoft = real(dzrf(Npulses,tb,pType,'ls',de,de));

if halfShift
  rfSoft = rfSoft.*((-1).^[0:Npulses-1]);
end

% using nominal peak b1, determine width of hard pulses
gambar = 4258; % Hz/g

% design the blip trapezoid
garea = 1/slsep/gambar*10; % mT/m*s, k-area of each blip
gzblip = dotrap(garea/10,gmax/10,gslew*100,dt).'*10; % mT/m

if ~minRFdur

    % matched-duration RF subpulses
    hpw = ceil(max(abs(rfSoft))./(2*pi*gambar*maxb1*dt));

    % interleave RF subpulses with gradient blips to form full pulses
    rf = kron(rfSoft(:),[ones(hpw,1);zeros(size(gzblip))]);
    NtSub = hpw+length(gzblip);
    
    g = repmat([zeros(hpw,1);gzblip(:)],[Npulses 1]);
    
    rf = rf./(sum(rf(:))*2*pi*gambar*dt)*sum(rfSoft); % convert to gauss

else

    % matched-amplitude RF subpulses
    rf = [];g = [];
    for ii = 1:Npulses
        hpw = ceil(abs(rfSoft(ii))./(2*pi*gambar*maxb1*dt));
        rftmp = rfSoft(ii)*ones(hpw,1);
        rftmp = rftmp./(sum(rftmp(:))*2*pi*gambar*dt)*sum(rfSoft(ii)); % convert to gauss
        rf = [rf; rftmp; zeros(size(gzblip))];
        g = [g; zeros(hpw,1); gzblip(:)];
    end
    
    NtSub = -1; 

end

switch pType
    case 'se'
        % remove the last blip, if 'se'
        g = g(1:end-length(gzblip)); % flip sign (gradient reversal to suppress fat) and remove last blip
        rf = rf(1:end-length(gzblip));
        ttipdown = dt*length(rf)/2;
    case {'st','ex'}
        % build a rewinder, if 'ex' or 'st'
        gzRew = dotrap(sum(g)*dt/2/10,gmax/10,gslew*100,dt).'*10; % mT/m
        g = g(1:end-length(gzblip)); % flip sign (gradient reversal to suppress fat) and remove last blip
        rf = rf(1:end-length(gzblip));
        ttipdown = dt*length(rf)/2;
        g = [g; -gzRew];
        rf = [rf; zeros(length(gzRew),1)];
end
