function res = mtimes(a,bb)

if streq(a.oneDorTwoD,'2D')
    if a.adjoint
        resAll = 0;
        Nc = size(a.sensMaps,2);
        for jj = 1:Nc
            % k -> image, via ifft
            kmat = zeros(size(a.kmask));
            kmat(a.kmask) = bb((jj-1)*length(bb)/Nc+1:jj*length(bb)/Nc);
            kmat = fftshift(kmat);
            
            if a.Nsh == 1 && max(abs(a.dphi)) == 0 && max(abs(a.dk)) == 0
                % no delays or phase shifts; just a single ifft2
                res = fftshift(ifft2(kmat))*(a.M*a.N);
            else
                fM2 = zeros(a.M,a.N/(2*a.Nsh),2*a.Nsh);
                for ii = 1:2*a.Nsh
                    if a.shotMask(mod(ii-1,a.Nsh)+1) % skip operations if shot was skipped.
                        % extract and DFT data for this segment
                        % go back to image domain; undo ifft2 weighting
                        fM2(:,:,ii) = ifft2(kmat(:,ii:2*a.Nsh:end))*(a.M*a.N/(2*a.Nsh));
                        % compensate k-space shift (since all ifft's assume that data starts at
                        % DC). We also apply k-space delays and phase shifts at this step
                        fM2(:,:,ii) = bsxfun(@times,fM2(:,:,ii),conj(a.pePhs(ii,:)));
                        fM2(:,:,ii) = bsxfun(@times,conj(a.fePhs(:,ii)),fM2(:,:,ii));
                    end
                end
                res = fftshift(reshape(col(reshape(fM2,[a.M*a.N/(2*a.Nsh) 2*a.Nsh])*a.C),[a.M a.N]));
            end
            
            if isfield(a,'immask')
                res = res(a.immask);
            else
                res = res(:);
            end
            resAll = resAll + res.*conj(a.sensMaps(:,jj));
        end
        res = resAll;
    else
        % image -> k, via fft
        resAll = [];
        for jj = 1:size(a.sensMaps,2) % loop over Rx coils
            
            tmp = zeros(numel(a.immask),1);
            tmp(a.immask) = bb.*a.sensMaps(:,jj);
            %bb = reshape(tmp,[a.M a.N]);
        
            if a.Nsh == 1 && max(abs(a.dphi)) == 0 && max(abs(a.dk)) == 0
                % no delays or phase shifts; just a single fft2
                res = fft2(fftshift(reshape(tmp,[a.M a.N])));
            else
                fM3 = reshape(reshape(fftshift(reshape(tmp,[a.M a.N])),[a.M*a.N/(2*a.Nsh) 2*a.Nsh])*a.C',...
                    [a.M a.N/(a.Nsh*2) a.Nsh*2]);
                % (apply conjugate delay and phase shifts here)
                res = zeros(a.M,a.N);
                for ii = 1:2*a.Nsh
                    if a.shotMask(mod(ii-1,a.Nsh)+1) % skip operations if shot was skipped.
                        % un-compensate k-space shifts
                        fM3(:,:,ii) = bsxfun(@times,a.fePhs(:,ii),fM3(:,:,ii));
                        fM3(:,:,ii) = bsxfun(@times,fM3(:,:,ii),a.pePhs(ii,:));
                        %fM3(:,:,ii) = fM3(:,:,ii).*(a.fePhs(:,ii)*a.pePhs(ii,:));
                        % go back to Fourier domain and
                        % interleave it with rest of k-space
                        res(:,ii:2*a.Nsh:end) = fft2(fM3(:,:,ii));
                    end
                end
            end
            
            res = fftshift(res);
            res = res(a.kmask);
            resAll = [resAll;res];
        end
        res = resAll;
    end
else % 1D
    if a.adjoint
        resAll = 0;
        Nc = size(a.sensMaps,2);
        for jj = 1:Nc
            % k -> image, via ifft
            kmat = zeros(size(a.kmask));
            kmat(a.kmask) = bb((jj-1)*length(bb)/Nc+1:jj*length(bb)/Nc);
            kmat = fftshift(kmat,2);
            
            if a.Nsh == 1 && max(abs(a.dphi)) == 0 && max(abs(a.dk)) == 0
                % no delays or phase shifts; just a single ifft
                res = ifftshift(ifft(kmat,[],2),2)*(a.N);
            else
                fM2 = zeros(a.M,a.N/(2*a.Nsh),2*a.Nsh);
                for ii = 1:2*a.Nsh
                    if a.shotMask(mod(ii-1,a.Nsh)+1) % skip operations if shot was skipped.
                        % extract and DFT data for this segment
                        % go back to image domain; undo ifft2 weighting
                        %                 fM2(:,:,ii) = ifft(kmat(:,ii:2*a.Nsh:end),[],2)*(a.M*a.N/(2*a.Nsh));
                        fM2(:,:,ii) = ifft(kmat(:,ii:2*a.Nsh:end),[],2)*(a.N/(2*a.Nsh));%*sqrt(a.M);
                        % compensate k-space shift (since all ifft's assume that data starts at
                        % DC). We also apply k-space delays and phase shifts at this step
                        %fM2(:,:,ii) = bsxfun(@times,fM2(:,:,ii),conj(a.pePhs(ii,:)));
                        %fM2(:,:,ii) = bsxfun(@times,conj(a.fePhs(:,ii)),fM2(:,:,ii));
                    end
                end
                fM2 = fM2.*conj(a.allPhs);
                res = ifftshift(reshape(col(reshape(fM2,[a.M*a.N/(2*a.Nsh) 2*a.Nsh])*a.C),[a.M a.N]),2);
            end
            
            if isfield(a,'immask')
                res = res(a.immask);
            else
                res = res(:);
            end
            resAll = resAll + res.*conj(a.sensMaps(:,jj));
        end
        
        res = resAll;
    else
        % image -> k, via fft
        resAll = [];
        
        for jj = 1:size(a.sensMaps,2) % loop over Rx coils
            
            tmp = zeros(numel(a.immask),1);
            tmp(a.immask) = bb.*a.sensMaps(:,jj);
            %bb = reshape(tmp,[a.M a.N]);
            
            if a.Nsh == 1 && max(abs(a.dphi)) == 0 && max(abs(a.dk)) == 0
                % no delays or phase shifts; just a single fft2
                res = fft(fftshift(reshape(tmp,[a.M a.N]),2),[],2);
            else
                fM3 = reshape(reshape(fftshift(reshape(tmp,[a.M a.N]),2),[a.M*a.N/(2*a.Nsh) 2*a.Nsh])*a.C',...
                    [a.M a.N/(a.Nsh*2) a.Nsh*2]);
                % (apply conjugate delay and phase shifts here)
                res = zeros(a.M,a.N);
                fM3 = fM3.*a.allPhs;
                for ii = 1:2*a.Nsh
                    if a.shotMask(mod(ii-1,a.Nsh)+1) % skip operations if shot was skipped.
                        % un-compensate k-space shifts
                        %fM3(:,:,ii) = bsxfun(@times,a.fePhs(:,ii),fM3(:,:,ii));
                        %fM3(:,:,ii) = bsxfun(@times,fM3(:,:,ii),a.pePhs(ii,:));
                        %fM3(:,:,ii) = fM3(:,:,ii).*(a.fePhs(:,ii)*a.pePhs(ii,:));
                        % go back to Fourier domain and
                        % interleave it with rest of k-space
                        res(:,ii:2*a.Nsh:end) = fft(fM3(:,:,ii),[],2);%*sqrt(a.M);
                    end
                end
            end
                        
            res = fftshift(res,2);
            res = res(a.kmask);
            resAll = [resAll;res];
        end
        res = resAll;
    end
end

