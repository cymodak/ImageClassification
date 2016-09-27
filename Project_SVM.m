close all
clear all
clc

%%
% Tunable Parameters
% numtrain: Number of trials for training
% numtest: Number of trials for testing
% Value of condition: Set i for 1-50, 50+i for 1-100, 100+i for 101-150...
%%
% Move to folder of code
if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

% Choose day
load 110215_1.mat
% Save structure in arrays and cells
condition = CDTTables{1,1}.condition;
start = CDTTables{1,1}.starttime;
stop = CDTTables{1,1}.stoptime;
times = CDTTables{1,1}.spikeTimes;
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

% Define bin-size
bin_size = 0.50;
spikearray = zeros(length(times),32,ceil(1/bin_size));

% Count number of spikes in each bin
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
spikearray1 = spikearray(:,:,floor((10*start_threshold)/(10*bin_size))+1:ceil(stop_threshold/bin_size));
spikearray = spikearray1;

% Choose day
% load 103015_1.mat
% % Save structure in arrays and cells
% condition = CDTTables{1,1}.condition;
% start = CDTTables{1,1}.starttime;
% stop = CDTTables{1,1}.stoptime;
% times = CDTTables{1,1}.spikeTimes;
% spikeTimes = zeros(length(times),32,399);
% electrode = CDTTables{1,1}.spikeElectrode;
% 
% for i = 1:1:length(times)
%     k = length(times{i,1});
%     for j = 1:1:k
%         last = length(times{i,1}{j,1});
%         spikeTimes(i,electrode{i,1}(j),1:last)  = times{i,1}{j,1};
%     end
% end
% 
% start_threshold = 0.3;
% stop_threshold = 0.65;
% spikeTimes(spikeTimes < start_threshold) = 0;
% spikeTimes(spikeTimes > stop_threshold) = 0;
% maxtime = max(max(max(spikeTimes)));
% 
% % for i = 1:1:length(times)
% %     for j = 1:1:32
% %         for k = 1:1:399
% %             if (spikeTimes(i,j,k) > stop(i))
% %                 spikeTimes(i,j,k) = spikeTimes(i,j,k)+1-stop(i);
% %             end
% %         end
% %     end
% % end
% 
% % Define bin-size
% bin_size = 0.50;
% spikearray = zeros(length(times),32,ceil(1/bin_size));
% 
% % Count number of spikes in each bin
% for i = 1:1:length(times)
%     for j = 1:1:32
%         for l=1:1:1/bin_size
%             [r1, c1] = find(spikeTimes(i,j,:)>(l-1)*bin_size);
%             [r2, c2] = find(spikeTimes(i,j,:)>(l)*bin_size);
%             spikearray(i,j,l) = length(r1) - length(r2);
%         end
%         %         for l=14:1:20
%         %             [r3, c3] = find(spikeTimes(i,j,:)>=1);
%         %             spikearray(i,j,12) = length(r3);
%         %         end
%     end
% end
% spikearray2 = spikearray(:,:,floor((10*start_threshold)/(10*bin_size))+1:ceil(stop_threshold/bin_size));

% spikearray = spikearray(:,:,2);
% 
% % Choose day
% load 1102.mat
% % Save structure in arrays and cells
% condition = CDTTables{1,1}.condition;
% start = CDTTables{1,1}.starttime;
% stop = CDTTables{1,1}.stoptime;
% times = CDTTables{1,1}.spikeTimes;
% spikeTimes = zeros(length(times),32,399);
% electrode = CDTTables{1,1}.spikeElectrode;
% 
% for i = 1:1:length(times)
%     k = length(times{i,1});
%     for j = 1:1:k
%         last = length(times{i,1}{j,1});
%         spikeTimes(i,electrode{i,1}(j),1:last)  = times{i,1}{j,1};
%     end
% end
% 
% start_threshold = 0.3;
% stop_threshold = 0.65;
% spikeTimes(spikeTimes < start_threshold) = 0;
% spikeTimes(spikeTimes > stop_threshold) = 0;
% maxtime = max(max(max(spikeTimes)));
% 
% % for i = 1:1:length(times)
% %     for j = 1:1:32
% %         for k = 1:1:399
% %             if (spikeTimes(i,j,k) > stop(i))
% %                 spikeTimes(i,j,k) = spikeTimes(i,j,k)+1-stop(i);
% %             end
% %         end
% %     end
% % end
% 
% % Define bin-size
% bin_size = 0.50;
% spikearray = zeros(length(times),32,ceil(1/bin_size));
% 
% % Count number of spikes in each bin
% for i = 1:1:length(times)
%     for j = 1:1:32
%         for l=1:1:1/bin_size
%             [r1, c1] = find(spikeTimes(i,j,:)>(l-1)*bin_size);
%             [r2, c2] = find(spikeTimes(i,j,:)>(l)*bin_size);
%             spikearray(i,j,l) = length(r1) - length(r2);
%         end
%         %         for l=14:1:20
%         %             [r3, c3] = find(spikeTimes(i,j,:)>=1);
%         %             spikearray(i,j,12) = length(r3);
%         %         end
%     end
% end
% spikearray3 = spikearray(:,:,floor((10*start_threshold)/(10*bin_size))+1:ceil(stop_threshold/bin_size));
% % spikearray = spikearray(:,:,2);
% 
% spikearray = cat(2,spikearray1,spikearray2);
% spikearray = spikearray(:,:,2);
% spikearray = spikearray1;
%Define number of training and testing samples
numtrain = 14;
numtest = 1;

% Build matrix of stimuli by clubbing all trials into an array

array = (1:15)';
accuracy = zeros(5,1);

count = 1;
while (count <= 5)
    psth1 = cell(50,1);
    psth_test1 = cell(50,1);
    for i = 1:1:50
        [r4, c4] = find (condition == 4*i-2);
        r4 = r4(array);
        for j = 1:1:numtrain
            temp = squeeze(spikearray(r4(j),:,:))';
            psth1{i} = [psth1{i} temp(:)];
        end
        for j = length(r4)-numtest+1:1:length(r4)
            temp = squeeze(spikearray(r4(j),:,:))';
            psth_test1{i} = [psth_test1{i} temp(:)];
        end
    end
    
    %% Training
    train = cell(1225,1);
    class = 1;
    
    for i = 1:1:49
        for j = i+1:1:50
            a = psth1{i};
            b = psth1{j};
            group = [ones(numtrain,1) ; zeros(numtrain,1)];
            c = [a' ; b'];
            train{class} = svmtrain(c,group,'kernel_function','mlp');
            class = class + 1;
        end
    end
    
    %% Testing
    % accuracy computes the name of correctly classified testing trials
    % accuracy = zeros(15,1);
    % count = 1;
    % Cell that stores the result of testing for each stimulus
    result_value = cell(50,1);
    for i = 1:1:50
        result = zeros(50,numtest);
        class = 1;
        for p = 1:1:49
            for q = p+1:1:50
                res = svmclassify(train{class},psth_test1{i}');
                for j = 1:1:numtest
                    if (res(j) > 0)
                        result(p,j) = result(p,j) + 1;
                    else
                        result(q,j) = result(q,j) + 1;
                    end
                end
                class = class + 1;
            end
        end
        temp = max(result);
        for s = 1:1:length(temp)
            [x1, y1] = find(result(:,s) == temp(s));
            [x2, y2] = find(x1 == i);
            accuracy(count) = accuracy(count) + length(x2);
        end
        result_value{i} = result;
    end
    array = circshift(array,10);
    count = count + 1;
end
            
            
