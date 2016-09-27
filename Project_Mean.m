close all
clear all
clc

% Move to folder of code
if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load 1019.mat
condition = CDTTables{1,1}.condition;
start = CDTTables{1,1}.starttime;
stop = CDTTables{1,1}.stoptime;
times = CDTTables{1,1}.spikeTimes;
unit = CDTTables{1,1}.spikeUnit;
spikeTimes = zeros(length(times),32,399);
electrode = CDTTables{1,1}.spikeElectrode;

for i = 1:1:length(times)
    k = length(times{i,1});
    for j = 1:1:k
        last = length(times{i,1}{j,1});
        spikeTimes(i,electrode{i,1}(j),1:last)  = times{i,1}{j,1};
    end
end

start_threshold = 0.3;
stop_threshold = 0.65;
spikeTimes(spikeTimes < start_threshold) = 0;
spikeTimes(spikeTimes > stop_threshold) = 0;
maxtime = max(max(max(spikeTimes)));

% for i = 1:1:length(times)
%     for j = 1:1:32
%         for k = 1:1:399
%             if (spikeTimes(i,j,k) > stop(i))
%                 spikeTimes(i,j,k) = spikeTimes(i,j,k)+1-stop(i);
%             end
%         end
%     end
% end

bin_size = 0.4;
spikearray = zeros(length(times),32,ceil(1/bin_size));

for i = 1:1:length(times)
    for j = 1:1:32
        for l=1:1:1/bin_size                 
            [r1, c1] = find(spikeTimes(i,j,:)>(l-1)*bin_size);
            [r2, c2] = find(spikeTimes(i,j,:)>(l)*bin_size);
            spikearray(i,j,l) = length(r1) - length(r2);
        end 
%         for l=14:1:20
%             [r3, c3] = find(spikeTimes(i,j,:)>=1);
%             spikearray(i,j,12) = length(r3);
%         end
    end
end
spikearray = spikearray(:,:,floor((10*start_threshold)/(10*bin_size))+1:ceil(stop_threshold/bin_size));
psth1 = cell(50,1);
psth_test1 = cell(50,1);
a = cell(50,1);

for i = 1:1:50
    [r4, c4] = find (condition == 4*i);
    psth1{i} = zeros(32,2);
    for j = 1:1:15
        temp = squeeze(spikearray(r4(j),:,:));
        psth1{i} = psth1{i} + temp;
    end
    a = psth1{i};
    a(:,1) = a(:,1)/(j*0.1);
    a(:,2) = a(:,2)/(j*0.25);
    psth1{i} = a;
%     psth1{i} = psth1{i}/(j*bin_size);
%     for j = 13
%         temp = squeeze(spikearray(r4(j),:,:));
%         psth_test1{i} = temp;
%     end
end
