function gr2ppe(fname,dur,ttipdown,g)

% fname: file name to write
% dur: duration in ms
% g: Ntx3 matrix of gradient sames, in mT/m (= 1/10 g/cm)

fp = fopen(fname,'wt'); % open file in text mode

fwrite(fp,sprintf('%f\n',dur)); % write duration in ms

fwrite(fp,sprintf('%f\n',max(abs(g(:,1)))));
fwrite(fp,sprintf('%f\n',max(abs(g(:,2)))));
fwrite(fp,sprintf('%f\n',max(abs(g(:,3)))));

fwrite(fp,sprintf('%f\n',ttipdown));

fwrite(fp,sprintf('%d\n',size(g,1)));

if max(abs(g(:,1))) > 0
  fwrite(fp,sprintf('%f\n',g(:,1)./max(abs(g(:,1)))));
else
  fwrite(fp,sprintf('%f\n',zeros(size(g,1),1)));
end
if max(abs(g(:,2))) > 0
  fwrite(fp,sprintf('%f\n',g(:,2)./max(abs(g(:,2)))));
else
  fwrite(fp,sprintf('%f\n',zeros(size(g,1),1)));
end
if max(abs(g(:,3))) > 0
  fwrite(fp,sprintf('%f\n',g(:,3)./max(abs(g(:,3)))));
else
  fwrite(fp,sprintf('%f\n',zeros(size(g,1),1)));
end

fclose(fp);



