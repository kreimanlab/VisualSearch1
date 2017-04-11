function oarr = myunique(iarr, varargin)

oarr=[];
for x = iarr
    if ~ismember(x, oarr)
        oarr = [oarr x];
    end
end
