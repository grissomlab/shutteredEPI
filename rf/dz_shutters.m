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
Nshots = 3; % number of shots/EPI segments
extraShotsForOverlap = 1; % # of extra shots to add so we get overlap.
cancelAlphaPhs = true;
% TODO: do we need to overencode if we add a shot? would need it to recon
% single shot well without overlap, but maybe not if we jointly recon them
%imFOV = 20.2; % cm, imaging FOV in shuttered dim
R = 4; % imaging acc factor
doSim = true; % do final Bloch simulation
if ~exist('inPlaneSimDim','var')
    inPlaneSimDim = [85 96];  % default in-plane sim dimensions. 96 is PE dim,
                              % which enables Nshots = 2,3,4 and R = 2,4
end
imFOV = 22;%0.2*inPlaneSimDim(2); % cm, imaging FOV in shuttered dim. 0.2 comes from res of B1+ maps
if ~exist('flip','var')
    flip = 65; % flip angle
end
flyback = false;
delayTolerance = 0; % fraction of nominal kslice width to pad so that
% the even or odd pulses (or both!) can be delayed without the RF falling
% down the ramp

tbw = [4 4]; % time bandwidth in slice, shutter dims
dthick = [0.4 imFOV/(R*Nshots)]; % slice thickness, shutter width (cm)

kw = tbw ./ dthick; % width of k-space coverage in each dimension (1/cm)

gz_area = kw(1) / 4257; % z(slice)-gradient area (g-s/cm)
gzmax = 2; % g/cm
gymax = 2; % g/cm
gslew = 20000; % g/cm/s
[gpos,ramppts] = domintrap(gz_area*(1+delayTolerance),gzmax,gslew,dt); % plateau sums to desired area
% remove last point since it is zero and will give two consecutive zeros in
% total waveform
gpos = gpos(1:end-1);
nFlyback = 0;
if flyback
    gzFlyback = dotrap(sum(gpos)*dt,gzmax,gslew,dt);
    gzFlyback = gzFlyback(1:end-1);
    gpos = [gpos -1*gzFlyback];
    nFlyback = length(gzFlyback);
end
Ntz = length(gpos);

% design slice-selective subpulse
rfSl = real(dzrf(round((Ntz-2*ramppts+1)/(1+delayTolerance))-nFlyback,tbw(1),'st','ls',0.01,0.01)); % arb units
% zero pad rf back to length of plateau if delayTolerance > 0
if delayTolerance > 0
    nPad = floor(((Ntz-2*ramppts+1)-length(rfSl))/2);
    rfSl = [zeros(1,nPad) rfSl zeros(1,nPad)];
    if length(rfSl) < Ntz-2*ramppts+1
        rfSl = [rfSl 0];
    end
end
% normalize to one radian flip
rfSl = rfSl./sum(rfSl);

% design the shutter envelope
if flip == 90
  if ~cancelAlphaPhs
    rfShut = real(dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'ex','ls',0.01,0.01)); % radians
  else
    [~,bShut] = dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'ex','ls',0.01,0.01);
    Bshut = ft(bShut);
    Bshut = Bshut.*exp(-1i*2*pi/round(kw(2)*Nshots*dthick(2))*1*(-round(kw(2)*Nshots*dthick(2))/2:round(kw(2)*Nshots*dthick(2))/2-1));
    bShut = ift(Bshut);
    aShut = b2a(bShut);
    bShut = ifft(fft(bShut).*exp(1i*angle(fft(aShut))));
    rfShut = real(b2rf(bShut));
  end
elseif flip == 180
  rfShut = real(dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'se','ls',0.01,0.01)); % radians
else % small-tip
  if ~cancelAlphaPhs
    rfShut = real(dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'st','ls',0.01,0.01)); % arb units
    % scale to target flip
    rfShut = rfShut./sum(rfShut)*flip*pi/180; % radians
  else
    bShut = dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'st','ls',0.01,0.01); % arb units
    Bshut = ft(bShut);
    Bshut = Bshut.*exp(-1i*2*pi/round(kw(2)*Nshots*dthick(2))*1*(-round(kw(2)*Nshots*dthick(2))/2:round(kw(2)*Nshots*dthick(2))/2-1));
    bShut = ift(Bshut);
    bShut = bShut*sind(flip/2);
    aShut = b2a(bShut);
    bShut = ifft(fft(bShut).*exp(1i*angle(fft(aShut))));
    rfShut = real(b2rf(bShut)); % radians
  end
end

%rfShut(1:2:end) = 0; % shut off every other pulse to image odd or even only

% construct the pulse with gaps for ramps
% flipping the rfSl for Even subpulses accommodates any off-centering of
% pulse due to earlier unequal zero padding
rfEPEven = kron(rfShut(2:2:end),...
    [zeros(1,2*ramppts+length(rfSl)-1+nFlyback) zeros(1,ramppts) fliplr(rfSl) zeros(1,ramppts-1+nFlyback)]);
rfEPOdd = kron(rfShut(1:2:end),...
    [zeros(1,ramppts) rfSl zeros(1,ramppts-1+nFlyback) zeros(1,2*ramppts+length(rfSl)-1+nFlyback)]);
if rem(length(rfShut),2)
    rfEPEven = [rfEPEven zeros(1,2*ramppts+length(rfSl)-1+nFlyback)];
    rfEPOdd = rfEPOdd(1:end-(2*ramppts+length(rfSl)-1+nFlyback));
end
rfEPEven = rfEPEven(1:end-nFlyback); % we will add half-area z rewinder later
rfEPOdd= rfEPOdd(1:end-nFlyback);
rfEP = rfEPEven + rfEPOdd;
ttipdown = length(rfEP)/2*dt*1000; % time into the pulse at which TE should start (ms) - calculate before we add rewinder zeros

% build total gz gradient waveform
if ~flyback
    gzEP = kron(ones(1,floor(length(rfShut)/2)),[gpos -gpos]);
    if rem(length(rfShut),2)
        gzEP = [gzEP gpos];
    end
else
    gzEP = repmat(gpos,[1 length(rfShut)]);
    gzEP = gzEP(1:end-nFlyback); % last rewinder will be half area
end

% get the gy blips
%gyBlip = dotrap((kw(2)/(length(rfShut)-1))/4257,gymax,gslew,dt);
gyBlip = dotrap(1/(Nshots*dthick(2))/4257,gymax,gslew,dt);
if rem(length(gyBlip),2)
    gyBlip = [gyBlip 0]; % add a zero to make it even length
end
if ~flyback
    % center gy blips between gz trapezoids
    % append zeros so that they straddle consecutive gz traps
    gyBlipPad = [zeros(1,Ntz-length(gyBlip)) gyBlip];
    gyEP = [zeros(1,length(gyBlip)/2) kron(ones(1,length(rfShut)-1),gyBlipPad)];
else
    % center gy blips on gz rewinders
    gyBlipPad = [zeros(1,Ntz-nFlyback+floor((nFlyback-length(gyBlip))/2)) ...
        gyBlip];
    gyBlipPad = [gyBlipPad zeros(1,Ntz-length(gyBlipPad))];
    gyEP = kron(ones(1,length(rfShut)-1),gyBlipPad);
end
gyEP = [gyEP zeros(1,length(gzEP)-length(gyEP))];

% calculate and add rewinders
gzRew = dotrap(sum(gpos(1:end-nFlyback))*dt/2,gzmax,gslew,dt);
if ~flyback
    gzEP = [gzEP ((-1)^rem(length(rfShut),2))*gzRew];
else
    gzEP = [gzEP -1*gzRew];
end
gyRew = -dotrap(sum(gyBlip)*dt*(length(rfShut)-1)/2,gymax,gslew,dt);
gyEP = [gyEP gyRew];

% zero pad waveforms to same length
gzEP = [gzEP zeros(1,max(length(gzEP),length(gyEP))-length(gzEP))];
gyEP = [gyEP zeros(1,max(length(gzEP),length(gyEP))-length(gyEP))];
gEP = [gyEP(:) gzEP(:)]; % stick them together into matrix
rfEP = [rfEP(:); zeros(size(gEP,1)-length(rfEP),1)];

% calculate FM waveform for slice-shifting
if ~flyback
    rfFM = repmat([ones(Ntz,1);-ones(Ntz,1)],[floor(length(rfShut)/2) 1]);
    if rem(length(rfShut),2)
        rfFM = [rfFM;ones(Ntz,1)];
    end
    rfFM = [rfFM; zeros(length(rfEP)-length(rfFM),1)];
else
    rfFM = ones(size(rfEP));
end

% calculate the phases to shift the slab to the other locations
phsMtx = angle(exp(1i*2*pi*(0:Nshots+extraShotsForOverlap-1)'/(Nshots+extraShotsForOverlap)*(0:length(rfShut)-1)));
rfPhs = kron(phsMtx,ones(1,Ntz));
rfPhs = rfPhs(:,1:end-nFlyback);
rfPhs = [rfPhs zeros(Nshots+extraShotsForOverlap,length(rfEP)-size(rfPhs,2))];

% confirm that the shutters go where we expect
figure;
y = (-64:63)/128*Nshots*dthick(2);
for ii = 1:Nshots+extraShotsForOverlap
    subplot(211)
    hold on
    plot(y,abs(fftshift(fft(kron(rfShut.*exp(1i*phsMtx(ii,:)),[1 0]),128))));
    ylabel 'approx flip angle'
    xlabel 'cm'
    subplot(212)
    hold on
    plot(y,sin(abs(fftshift(fft(kron(rfShut.*exp(1i*phsMtx(ii,:)),[1 0]),128)))));
    ylabel 'approx |Mxy|'
    xlabel 'cm'
end
title 'All shutters'

flips = 0:5:90;
TR = 500; % ms; this is really the time between shots
T1 = 2000; % ms 7T grey matter
nReps = 10; % number of TRs to run
signal = zeros(length(flips), 1);
for ll = 1 : length(flips)
    flipAngles = zeros(length(y), Nshots+extraShotsForOverlap);
    for ii = 1 : Nshots+extraShotsForOverlap
        flipAngles(:, ii) = flips(ll) / 90 * abs(fftshift(fft(kron(rfShut.*exp(1i*phsMtx(ii,:)),[1 0]),128)));
    end
    Mz = ones(length(y), 1);
    Mxy = zeros(length(y), Nshots+extraShotsForOverlap);
    for ii = 1 : nReps
        for jj = 1 : Nshots+extraShotsForOverlap
            % apply pulse
            Mxy(:, jj) = Mz .* sin(flipAngles(:, jj));
            Mz = Mz .* cos(flipAngles(:, jj));
            % relax
            Mz = 1 - (1 - Mz) * exp(-TR / T1);
        end
%         figure;
%         subplot(211)
%         plot(y, abs(Mxy));
%         axis([y(1) y(end) 0 1]);
%         subplot(212)
%         plot(y, Mz);
%         axis([y(1) y(end) -1 1]);
%         drawnow;
%         pause;
    end
    signal(ll) = mean(sqrt(sum(Mxy.^2, 2)), 1);
end
figure;plot(flips,signal);

signalsFullEx = zeros(length(flips), 1);
for ll = 1 : length(flips)
    flipAngles = zeros(length(y), Nshots+extraShotsForOverlap);
    for ii = 1 : Nshots+extraShotsForOverlap
        flipAngles(:, ii) = flips(ll) * pi / 180;
    end
    Mz = ones(length(y), 1);
    Mxy = zeros(length(y), Nshots+extraShotsForOverlap);
    for ii = 1 : nReps
        for jj = 1 : Nshots+extraShotsForOverlap
            % apply pulse
            Mxy(:, jj) = Mz .* sin(flipAngles(:, jj));
            Mz = Mz .* cos(flipAngles(:, jj));
            % relax
            Mz = 1 - (1 - Mz) * exp(-TR / T1);
        end
    end
    signalsFullEx(ll) = mean(sqrt(sum(Mxy.^2, 2)), 1);
end
hold on
plot(flips, signalsFullEx);
xlabel 'degrees'
ylabel 'Mean ssq\{|Mxy|\}'
legend('Shuttered', 'FullEX')
axis([1 90 0 1])
grid on 

% plot the pulses
figure
subplot(411)
t = (0:length(rfEP)-1)*dt*1000; % time in ms
plot(t,rfEP);
xlabel 'ms',ylabel 'radians'
title 'RF pulse'
subplot(412)
plot(t,gEP);
xlabel 'ms',ylabel 'g/cm'
title 'Gradient waveforms'
legend('gz','gy');
subplot(413)
plot(t,rfFM)
xlabel 'ms',ylabel 'arb units'
title 'Normalized FM waveform for slice shifting'
subplot(414)
plot(t,rfPhs)
xlabel 'ms',ylabel 'Radians'
title 'Phase waveforms to shift shutter to each location (one per shot)'
c = axis; axis([c(1) c(2) -4 4]);

if doSim
    disp 'Bloch-simulating final pulses'
    [mxy,~,alpha,beta] = blochsim_spinor(rfEP/(2*pi*4257*dt),gEP,...
        [dthick(2)*Nshots 10*dthick(1)],[128 128],zeros(128),dt);
    mxy = mxy.';
    beta = beta.';
    if flip <= 90
        figure;im((-64:63)/128*10*dthick(1),(-64:63)/128*dthick(2)*Nshots,mxy);
        ylabel 'y (shutter), cm'
        xlabel 'z (slice), cm'
        title 'Excitation profile'
    else
        figure;im((-64:63)/128*10*dthick(1),(-64:63)/128*dthick(2)*Nshots,abs(beta).^2);
        ylabel 'y (shutter), cm'
        xlabel 'z (slice), cm'
        title 'Refocusing profile'
    end
    mxyInPlane = zeros([inPlaneSimDim Nshots+extraShotsForOverlap]); % to match XY's sims
    for ii = 1:Nshots+extraShotsForOverlap
        mxyInPlane(:,:,ii) = blochsim_spinor(rfEP.*exp(1i*rfPhs(ii,:)')/(2*pi*4257*dt),[0*gEP(:,2) gEP(:,1)],...
            [20.2 20.2],inPlaneSimDim,zeros(inPlaneSimDim),dt);
    end
    save(sprintf('dz_shutters_tb%d_nshots%d_R%d',tbw(2),Nshots,R),'mxyInPlane');
end

% save pulses for LT pTx design testing
save(sprintf('dz_shutters_tb%d_nshots%d_R%d_pulses',tbw(2),Nshots,R),'rfEP','gEP','dt','dthick','tbw');

% write to philips PPE
rfAng = sum(rfEP)*180/pi; % convert to degrees
rf2ppe('shutter_rf.txt',length(rfEP)*dt*1000,rfAng,ttipdown,0,rfEP,rfFM);
gr2ppe('shutter_grads.txt',size(gEP,1)*dt*1000,ttipdown,10*[zeros(size(gEP,1),1) gEP]);
