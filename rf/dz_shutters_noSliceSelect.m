% Shutter EPI excitations.
% Assumes 90 degree tip;
% lower flips should switch to 'st' design of rfShut,
% and scale flip accordingly.

% total # of shutters is R*Nshots
% so the width of each shutter should be FOV/(R*Nshots),
% and the phase encode XFOV should be (Nshots)*(shutter width)

% rfSl: The slice-selective pulse (radians; scaled to 1)
% rfShut: The shutter envelope (radians; scaled to pi/2)
% gpos: The slice-select trapezoid (g/cm)
% gyBlip: The blip between the subpulses (g/cm)
% rfEP: The total constructed pulse (radians)
% gzRew: the slice-select rewinder
% gyRew: the phase encode/slab-select rewinder
% ttipdown: the time into the pulse when the TE period should start
% rfFM: normalized FM waveform for slice shifting
% rfPhs: matrix of RF phase waveforms that shift the shutter for each shot
% phsMtx: matrix of RF phases (one per subpulse) that shift the shutter for
%         each shot

% dependencies:
% - JP's RF toolbox (for dzrf)
% - domintrap.m - design trapezoid for a target area under the flat
%   portion (for non-ramp-sampled pulses).
% - dotrap.m - design trapezoid for total target area (for rewinders)
% - Fessler IRT for im (could replace with imagesc)
% - rf2ppe.m - write rf pulse to Philips PPE
% - gr2ppe.m - write gradient pulses to Philips PPE

dt = 6.4e-6; % s, dwell times of RF + gradient pulses
Nshots = 4; % number of shots/EPI segments
%imFOV = 20.2; % cm, imaging FOV in shuttered dim
R = 3; % imaging acc factor
doSim = true; % do final Bloch simulation
maxB1 = 0.1; % gauss, max b1 amplitude
gymax = 4; % g/cm
gslew = 20000; % g/cm/s
if ~exist('inPlaneSimDim','var')
    inPlaneSimDim = [85 96];  % default in-plane sim dimensions. 96 is PE dim,
                              % which enables Nshots = 2,3,4 and R = 2,4
end
imFOV = 0.2*inPlaneSimDim(2); % cm, imaging FOV in shuttered dim
if ~exist('flip','var')
    flip = 90; % flip angle
end
if flip < 90
    pType = 'st';
else
    pType = 'ex';
end
flyback = false;
delayTolerance = 0; % fraction of nominal kslice width to pad so that
% the even or odd pulses (or both!) can be delayed without the RF falling
% down the ramp

tbw = [4 4]; % time bandwidth in slice, shutter dims
dthick = [0.1 imFOV/(R*Nshots)]; % slice thickness, shutter width (cm)

kw = tbw ./ dthick; % width of k-space coverage in each dimension (1/cm)

% use PINS pulse design code to design the waveforms; move the PINS
% z-gradient to the y-axis
[rfEP,gEP,ttipdown,NSub,NtSub] = dz_pins(tbw(2),Nshots*dthick(2),dthick(2)*10,0.01,maxB1,gymax*10,gslew/100,...
    dt,false,false,pType);
% NSub = # of subpulses
% NtSub = # of time points in a subpulse + gy blip pair
gEP = gEP/10; % back to gauss/cm
ttipdown = ttipdown*1000; % ms

% calculate the phases to shift the slab to the other locations
phsMtx = angle(exp(1i*2*pi*(0:Nshots-1)'/Nshots*(0:NSub-1)));
rfPhs = kron(phsMtx,ones(1,NtSub));
rfPhs = [rfPhs zeros(Nshots,length(rfEP)-size(rfPhs,2))].';

% plot the pulses
figure
subplot(311)
t = (0:length(rfEP)-1)*dt*1000; % time in ms
plot(t,rfEP);
xlabel 'ms',ylabel 'radians'
title 'RF pulse'
subplot(312)
plot(t,gEP);
xlabel 'ms',ylabel 'g/cm'
title 'Gradient waveforms'
legend('gy');
subplot(313)
plot(t,rfPhs)
xlabel 'ms',ylabel 'Radians'
title 'Phase waveforms to shift shutter to each location (one per shot)'
c = axis; axis([c(1) c(2) -4 4]);

% write to philips PPE - the first one is real-valued
rfAng = sum(rfEP)*2*pi*4258*dt*180/pi; % convert to degrees
rf2ppe('shutter_rf_shot1.txt',length(rfEP)*dt*1000,rfAng,ttipdown,0,rfEP);
for ii = 2:Nshots
    rfShift = rfEP.*exp(1i*rfPhs(:,ii));
    rfFM = diff([0;unwrap(angle(rfShift))])/dt/2/pi;
    rfAng = sum(abs(rfShift))*2*pi*4258*dt*180/pi; % convert to degrees
    rf2ppe(sprintf('shutter_rf_shot%d.txt',ii),length(rfEP)*dt*1000,rfAng,ttipdown,0,abs(rfEP),rfFM);
end
gr2ppe('shutter_grads.txt',size(gEP,1)*dt*1000,ttipdown,10*[zeros(size(gEP,1),1) gEP zeros(size(gEP,1),1)]);
