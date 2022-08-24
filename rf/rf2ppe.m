function rf2ppe(fname,dur,flipangle,ttipdown,phs,rf_am,rf_fm)

% fname: string for filename
% dur: duration in ms
% ttipdown: time to tip-down
% flipangle: degrees
% rf: rf pulse
% phs: initial phase of pulse

fp = fopen(fname,'wt'); % open file in text mode

fwrite(fp,sprintf('%f\n',dur)); % write duration

fwrite(fp,sprintf('%f\n',flipangle)); % write flip angle (= integral of rf_am*gam*2*pi*dt)

fwrite(fp,sprintf('%f\n',phs)); % write initial phase

if exist('rf_fm');
  fwrite(fp,sprintf('%f\n',max(abs(rf_fm)))); % write the peak fm mod in Hz
else
  fwrite(fp,sprintf('%f\n',0));
end

fwrite(fp,sprintf('%f\n',ttipdown));

fwrite(fp,sprintf('%d\n',length(rf_am))); % write number of samples

% let the am samples go from -1 to 1
if max(abs(rf_am)) > 0
  fwrite(fp,sprintf('%f\n',rf_am./max(abs(rf_am))));
else
  fwrite(fp,sprintf('%f\n',zeros(length(rf_am),1)));
end

if nargin == 6
  % write zeros for fm modulation
  fwrite(fp,sprintf('%f\n',zeros(length(rf_am),1)));
else
  fwrite(fp,sprintf('%f\n',rf_fm./max(abs(rf_fm))));
end

fclose(fp);
