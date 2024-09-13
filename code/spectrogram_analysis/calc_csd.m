% % Calculate the CSD
% 
% inputs
% - LFP (channels x time x trials)
% - h (which neighbor)
% - spat_fs (spatial sampling)
% 
% output
% - CSD: cell 1 x trials, with arrays of size channels x time (fieldtrip format)

function CSD = calc_csd(LFP,h,spat_fs)

    CSD= cell(1,size(LFP,3));
    for t = 1:size(LFP,3)
        CSD{t} = zeros(size(LFP,1)-2*h,size(LFP,2));
    
        for c = 1+h:size(LFP,1) - h
            CSD{t}(c-h,:) = (LFP(c-h,:,t) - LFP(c,:,t) * 2 + LFP(c+h,:,t))/(h^2 * (1/spat_fs)^2) ;
        end
    end

end