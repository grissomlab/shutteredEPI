function [s,x] = epiPhsDelayEst_nd3(yref,ytest,stopThresh,delayRange)

if nargin < 3
    stopThresh = 0.00001;
end

[N,L,C] = size(yref);   % N: # spatial locs/k-space points; C: # coils
                        % L: # of k-space lines to align with the reference

% stage 0: find a coarse delay
delays = -delayRange/2:delayRange/2;
for ii = 1:length(delays)
    corr(ii) = sum(col(abs(yref).*circshift(abs(ytest),delays(ii))));
end
[~,ind] = max(corr);
% convert the best shift into a slope
s = delays(ind);
%ytest = circshift(ytest,delays(ind));
% get initial phase difference as mean phase of inner product - REDFLAG may
% have wrong sign
s(2) = -angle(sum(col(conj(yref).*circshift(ytest,delays(ind)))));s = s(:);

% FT both signals to image domain
yref = fftshift(ifft(fftshift(yref,1),[],1),1);
yref = yref(:);
ytest = fftshift(ifft(fftshift(ytest,1),[],1),1);
ytest = ytest(:);
% phase coordinates
x = 2*pi*((0:N-1)'/N - 1/2);
x = repmat(x,[1 L C]);x = x(:); % duplicate to all test lines
% signal magnitude product and phase difference - calculate a priori
t1 = abs(yref.*ytest);
sigPhsDiff = angle(ytest.*conj(yref));

% stage II: update both slope and DC phase shift

costOld = Inf;
%s = [0;s];
A = [x 0*x+1];

costNew = sum(abs(yref-ytest.*exp(1i*(A*s))).^2);

while costNew < (1-stopThresh)*costOld
    
    costOld = costNew;
    
    % update s
    t2 = A*s + sigPhsDiff;
    t3 = sin(t2);
    t4 = mod(t2+pi,2*pi)-pi;
    t5 = t1.*t3;
    s = s - (A' * bsxfun(@times, t5 ./ (t4 + eps), A) + 0.00001 * eye(size(A, 2))) \ (A' * t5);
    
    % recalculate cost
    costNew = sum(abs(yref-ytest.*exp(1i*(A*s))).^2);
    
end


