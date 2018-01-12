function [ SigPars ] = RandomMixingSignal(MetaPars)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

rootFreq_inds = sort(randsample(length(MetaPars.f)/2,2),'descend');
rootFreqs = MetaPars.f(rootFreq_inds);
SigPars = [ rootFreqs(1)    ... 1. high root frequency
            rootFreqs(2)    ... 2. low root frequency
            [rand(1) + 0.25]       ... 3. high root amplitude
            [rand(1) + 0.25]        ... 4. low root amplitude
            [rand(1) + 0.5]         ... 5. A1
            [rand(1) + 0.2]         ... 6. A2
            [rand(1) + 0.1]         ... 7. A3
            rand*2         ... 8. high root half-bandwidth
            rand*2         ... 9. low root half-bandwidth
            ceil(rand*3)    ... 10. high root n components
            ceil(rand*3)    ... 11. low root n components
            rand*0.5          ... 12. noise scale factor
            round([rand(1,1)*50+2])    ...13. added (non-mixing) component 1 freq
            [rand(1)/2 + 0.05]                    ...14. added (non-mixing) component 1 amp
            rand*1                    ...15. added (non-mixing) component 1 half-bandwidth
            round([rand(1,1)*50+2])    ...16. added (non-mixing) component 2 freq
            [rand(1)/2 + 0.05]                    ...17. added (non-mixing) component 2 amp
            rand*2                    ...18. added (non-mixing) component 2 half-bandwidth
            round([rand(1,1)*50+2])    ...19. added (non-mixing) component 3 freq
            [rand(1)/2 + 0.05]                   ...20. added (non-mixing) component 3 amp
            rand*4];%                    %21. added (non-mixing) component 3 half-bandwidth
end

