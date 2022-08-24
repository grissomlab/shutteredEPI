function sz = size(a,dim)

sz = [sum(a.kmask(:))*size(a.sensMaps,2) numel(a.kmask)];

if exist('dim','var')
    sz = sz(dim);
end
