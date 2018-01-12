function [obj] = FM_FFT_S_full(obj)

N_freqs = obj.N_freqs;
phase_bins = linspace(-pi,pi,6);
data  = angle(obj.data);
data1 = angle(obj.data);
data2 = angle(obj.data);
data3 = angle(obj.data);

% this is to randomize phase on one of three indices
for i=1:size(data,2) %for each frequency
    %k=randsample(3,1);
    %if k==1 
        data1(:,i) = data1(randsample(size(data,1),size(data,1)),i);
    %end
%     if k==2
%         data2(:,i) = data2(randsample(size(data,1),size(data,1)),i);
%     end
%     if k==3 
%         data3(:,i) = data3(randsample(size(data,1),size(data,1)),i);
%     end
%     
%     elseif k==2
%         data2(:,i) = data1(randsample(size(data,1),size(data,1)),i);
%     else
%         data3(:,i) = data1(randsample(size(data,1),size(data,1)),i);
%     end
end
        
obj.fft_pII_analysis.S_full = NaN(N_freqs,N_freqs,N_freqs,11);
FFT_inds = [4 2 2;... #1:  1 -1 -1
            4 2 1;... #2:  1 -1 -2
            4 1 4;... #3:  1 -2  1
            4 1 2;... #4:  1 -2 -1
            5 2 1;... #5:  2 -1 -2
            5 1 2;... #6:  2 -2 -1
            4 1 1;... #7:  1 -2 -2
            4 1 5];  %#8:  1 -2  2

for i=1:N_freqs
    i
    for j = 1:i
        for k = 1:j
            temp_binned_data1 = histcn([data1(:,i),data2(:,j),data3(:,k)],...
                phase_bins, phase_bins, phase_bins);
            temp_binned_data2 = (temp_binned_data1 - mean(temp_binned_data1(:)))/sum(temp_binned_data1(:));
        
            S__ = fftn(temp_binned_data2,[5 5 5]);
            S_ = fftshift(S__);
            S = abs(S_).^2;
            for w=1:8
                obj.fft_pII_analysis.S_full(i,j,k,w) = S(FFT_inds(w,1),FFT_inds(w,2),FFT_inds(w,3));
            end
            temp1 = S(3,:,:); obj.fft_pII_analysis.S_full(i,j,k,9) = max(temp1(:));
            temp2 = S(:,3,:); obj.fft_pII_analysis.S_full(i,j,k,10) = max(temp2(:));
            temp3 = S(:,:,3); obj.fft_pII_analysis.S_full(i,j,k,11) = max(temp3(:));
        end
    end
end

end