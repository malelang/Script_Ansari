function [pks,pk_locs,tr_locs] = myfindpeaks(X,varargin)

    [pks,pk_locs] = findpeaks(X,varargin{:});
    
    if(nargout > 2)
        tr_locs = zeros(length(pk_locs)-1,1);
        for i = 1:length(pk_locs)-1
            [~,mi] = min(X(pk_locs(i)+1:pk_locs(i+1)));
            tr_locs(i) = pk_locs(i) + mi;
        end    
    end

end

