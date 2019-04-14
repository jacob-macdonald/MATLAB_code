% This function characterizes the gating across a PC-VIPR scan. Missed
% beats are counted and corrected for. The variation in heart rate across
% the corrected ECG is then plotted.

% Author: Jacob Macdonald
% Date: October 3, 2017

clear all
close all

%% Open gating track

% ******EDIT THIS VARIABLE******
name = 'Gating_Track_y2017_m4_d27_h16_m14_s43.pcvipr_track';

disp(['Opening - ',name])
if name(end-4:end)=='track'
    
    fid = fopen(name);
    gate = fread(fid,'int32','b');
    gate = reshape(gate,[numel(gate)/4 4]);
    gate(:,5) = 1:size(gate,1);
    gate = sortrows(gate,3);
    fclose(fid);
    
    
    data.ecg = gate(:,1);
    data.resp= 4095-gate(:,2);
    data.time =gate(:,3)/1e6;
    data.prep= gate(:,4);
    data.acq = gate(:,5);
    
else
    
    fid = fopen(name);
    raw = fread(fid,'int32','b');
    raw = reshape(raw,[numel(raw)/5 5]);
    fclose(fid);
    
    data.ecg = raw(:,1);
    data.resp= 4095-raw(:,2);
    data.time =raw(:,3)/1e6;
    data.prep= raw(:,4);
    data.acq = raw(:,5);
    
end

%% Assess original gating

N = length(data.ecg);
count = 1;

for n = 2:1:N
    if(data.ecg(n) < data.ecg(n-1))
        RR_old(count) = data.ecg(n-1);
        count = count + 1;
    end
end

% Note this median value will likely disagree with that reported by the
% PC-VIPR recon, as that filters out projections with ECG>5000ms in advance
% of its calculation. We don't do that here, as those projections may
% result from multiple missed heart beats.

M1 = median(RR_old);
HB_old = round(1/(M1/60000));

disp(' ')
disp('Original gating statistics:')
Text = sprintf('    Median RR interval: %d ms.',M1);
disp(Text)
Text = sprintf('    Median heart rate: %d bpm.',HB_old);
disp(Text)

%% Clean gating signal

diff = zeros(N-1,1);

for n=1:N-1
    diff(n) = data.ecg(n+1) - data.ecg(n);
end

% Will fail in highly corrputed gating with no notable trends
readout_length = median(diff);

% This loop removes noise spikes in the data
for n = 1:N-1
    if (diff(n)>(readout_length+25)) || (diff(n)<(readout_length-25)) % make this more robust later
        if diff(n)<0
            data.ecg(n+1) = data.ecg(n+1); % if diff is because of end of RR interval, ignore
        else
            data.ecg(n+1) = round(mean([data.ecg(n) data.ecg(n+2)])); % avg for increased robustness in the presence of variable RR lengths
            if(data.ecg(n+1)<data.ecg(n)) % this deals with any spikes at the end of an RR interval
                data.ecg(n+1)=data.ecg(n)+readout_length;
            end
        end
    end
end

%******* EDITING NOTES: Script seems good up to here. Spike detection could
%still be made more robust on a case-by-case basis, but not the concern
%now. ********

% Recheck the RR values now that extreme outliers have been removed
count = 1;
peaks = zeros(N,1);
for n = 2:1:N
    if(data.ecg(n) < data.ecg(n-1))
        RR_mid(count) = data.ecg(n-1);
        peaks(n-1) = 1;
        count = count + 1;
    end
end

%typical_RR = median(RR_mid);
index = zeros(count,1);
index(2:count) = find(peaks);
L = length(index);
scale_vector = zeros(L-1,1);

str = 'j'; % random string value for first loop
temp = data.ecg; % temp holder for data.ecg in case reset is needed for second loop

while str ~='y';
    for l = 2:L
        if l<25 %l<50
            %RR_window = linspace(1,100,100);
            RR_window = linspace(1,50,50);
        elseif l>L-26 %l > L-51
            %RR_window = linspace(L-100,L-1,100);
            RR_window = linspace(L-50,L-1,50);
        else
            %RR_window = linspace(l-49,l+50,100);
            RR_window = linspace(l-24,l+25,50);
        end
        
        sliding_RR = median(RR_mid(RR_window));
        
        if str == 'n'
            sliding_RR = RR_hardcode;
        end
        
        scale = data.ecg(index(l)) / sliding_RR;
        scale_vector(l-1) = round(scale);
        
        if scale < 0.75 && l < L % account for early false triggers
            if (data.ecg(index(l))+data.ecg(index(l+1))) < 1.25*sliding_RR
                data.ecg((index(l)+1):index(l+1)) = data.ecg((index(l)+1):index(l+1)) + data.ecg(index(l));
            end
        elseif scale_vector(l-1) > 1 % account for missed heart beats
            range = index(l) - index(l-1);
            subrange = ceil(range/scale_vector(l-1));
            last_subrange = range - (scale_vector(l-1)-1)*subrange;
            
            for n = 1:scale_vector(l-1)
                if n == scale_vector(l-1)
                    data.ecg((index(l-1)+1+(n-1)*subrange):(index(l-1)+range)) = data.ecg((index(l-1)+1):(index(l-1)+last_subrange));
                else
                    data.ecg((index(l-1)+1+(n-1)*subrange):(index(l-1)+n*subrange)) = data.ecg((index(l-1)+1):(index(l-1)+subrange));
                end
            end
        end
    end


    % Recheck the RR values now that missed heart beats have been adjusted
    count = 1;
    peaks = zeros(N,1);
    data_norm = zeros(N,1);
    temp_index = 1;
    for n = 2:1:N
        if(data.ecg(n) < data.ecg(n-1))
            RR_new(count) = data.ecg(n-1);
            time(count) = data.time(n-1);
            peaks(n-1) = 1;
            count = count + 1;
            data_norm(temp_index:(n-1)) = data.ecg(temp_index:(n-1)) / data.ecg(n-1);
            temp_index = n;
        end
    end
   
    
    % Plot mean +std dev for cleaned RR data
    
    M = length(RR_new);
    for ii = 1:M
        if ii<25 %ii<50
            %RR_window = linspace(1,100,100);
            RR_window = linspace(1,50,50);
        elseif ii>M-26 %l > L-51
            %RR_window = linspace(L-100,L-1,100);
            RR_window = linspace(M-50,M-1,50);
        else
            %RR_window = linspace(l-49,l+50,100);
            RR_window = linspace(ii-24,ii+25,50);
        end
        
        RR_vector(ii) = median(RR_new(RR_window));
        std_vector(ii) = std(RR_new(RR_window));
    end
    
        data_norm = data_norm * median(RR_vector);
    
    figure()
    x = 1:M;
    fill([x fliplr(x)],[std_vector+RR_vector fliplr(RR_vector-std_vector)],[.9 .9 .9],'linestyle','none')
    ylim([0 2000])
    ylabel('RR interval length [ms]')
    title('RR time response')
    hold on
    plot(RR_vector, 'LineWidth',2)
    plot(RR_new)
    hold off
    
    M2 = median(RR_new);
    HB_new = round(1/(M2/60000));
    
    disp(' ')
    disp('Cleaned gating statistics:')
    Text = sprintf('    Median RR interval: %d ms.',M2);
    disp(Text)
    Text = sprintf('    Median heart rate: %d bpm.',HB_new);
    disp(Text)
    
    % Characterize missed heart beats
    total_hb = sum(scale_vector);
    missed_hb = scale_vector(find(scale_vector>1))-1; % vector of how many successive HB were missed each miss
    total_missed_hb = sum(missed_hb);
    percent_missed = total_missed_hb/total_hb*100;
    consecutive_misses = sum(missed_hb-1);
    
    disp(' ')
    disp('Missed heart beat statistics:')
    Text = sprintf('    Total heart beats: %d',total_hb);
    disp(Text);
    Text = sprintf('    Missed heart beats: %d (%2.1f%%)',total_missed_hb, percent_missed);
    disp(Text);
    Text = sprintf('    Consecutively missed heart beats: %d',consecutive_misses);
    disp(Text);
    disp(' ')
    
    
    % The code below plots before and after RR histograms and saves the cleaned
    % gating file.
    
    figure()
    subplot(1,2,1)
    histogram(RR_old, [0:10:3000], 'Normalization', 'probability','EdgeColor','none')
    title('Original RR distribution')
    xlabel('RR interval [ms]')
    ylabel('Probability')
    subplot(1,2,2)
    histogram(RR_new, [0:10:3000], 'Normalization', 'probability','EdgeColor','none')
    title('Cleaned RR distribution')
    xlabel('RR interval [ms]')
    ylabel('Probability')

    str = input('Are you satisfied with the gating? (y/n) ', 's');
    if str == 'n'
        data.ecg = temp;
        RR_hardcode = input('Input RR length for new cleaning: ');
    end
end

str = input('Save new gating file? (y/n) ', 's');
if str == 'y'

    str = input('Variable (v) or normalized (n) gating? (v/n) ', 's');
    if str == 'v'
        %Export cleaned gating in new file
        name2 = strcat('Clean_',name);
        gate2 = gate;
        gate2(:,1) = data.ecg(:);
        
        fid = fopen(name2,'w');
        gate2 = sortrows(gate2,5);
        gate2 = gate2(:,1:4);
        gate2 = reshape(gate2, [numel(gate2) 1]);
        fwrite(fid,gate2,'int32','b');
        fclose(fid);
    elseif str == 'n'
        name2 = strcat('Normalized_',name);
        gate2 = gate;
        gate2(:,1) = data_norm(:);
        
        fid = fopen(name2,'w');
        gate2 = sortrows(gate2,5);
        gate2 = gate2(:,1:4);
        gate2 = reshape(gate2, [numel(gate2) 1]);
        fwrite(fid,gate2,'int32','b');
        fclose(fid);
    end
end