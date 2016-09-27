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

% Load Particular Day
load 1019.mat
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

% Windowing
start_threshold = 0.2;
stop_threshold = 0.75;
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
bin_size = 0.05;
spikearray = zeros(length(times),32,ceil(1/bin_size));

% Count number of spikes in each bin
for i = 1:1:length(times)
    for j = 1:1:32
        for l=1:1:ceil(1/bin_size)
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
spikearray = spikearray(:,:,floor(start_threshold/bin_size)+1:ceil(stop_threshold/bin_size));
psth1 = cell(50,1);
psth_test1 = cell(50,1);

%Define number of training and testing samples
numtrain = 14;
numtest = 1;

% Build matrix of stimuli by clubbing all trials into an array
for i = 1:1:50
    [r4, c4] = find (condition == i);
    for j = 1:1:numtrain
        temp = squeeze(spikearray(r4(j),:,:))';
        psth1{i} = [psth1{i} temp(:)];
    end
    for j = length(r4)-numtest+1:1:length(r4)
        temp = squeeze(spikearray(r4(j),:,:))';
        psth_test1{i} = [psth_test1{i} temp(:)];
    end
end

% Store thresholds and locations of split for each strong classifier in
% corresponding cells
% Each cell contains locations of split for each weak classifier
threshold_value = cell(1225,1);
test_value = cell(1225,1);
sign_final_value = cell(1125,1);

%%
class = 1;
feig = 4;
for p = 1:1:49
    for q = p+1:1:50
        train_data = [psth1{p} psth1{q}];
        % for i = 1:1:24
        %     train_data(:,i)= train_data(:,i) - mean(train_data,2);
        % end
        [U, S, V] = svds(train_data,feig);
        w1 = U'*train_data;        
        wtrain = w1';
        
        % imvalue: scores the sorted weights for each eigenvector
        % imindex: book-keeping matrix that maps sorted  weights to original weights
        % Probability: Updated at the end of each iteration
        % alpha: Weight of each iteration
        % h: Classification label at the end of each iteration
        % sign_final: Direction of classification for each iteration
        % totalerror: A metric for measuring training accuracy
        % Sorting weights in a columnwise manner
        [imvalue,imindex] = sort(wtrain,1);        
        threshold = zeros(size(imindex));
        
        % Gives time 't'
        iter = 1;
        % Signals the end of time period
        maxiter = 20;
        
        prob = 1/size(imindex,1).*ones(size(imindex,1),maxiter);
        alpha = zeros(maxiter,1);
        h = ones(size(wtrain,1),maxiter);
        y = [ones(numtrain,1); -1*ones(numtrain,1)];
        test = zeros(maxiter,2);
        totalerror = 0;
        sign_final = zeros(maxiter,1);
        
        while(iter <= maxiter)
            error1 = zeros(size(imindex));
            sign = zeros(size(imindex));
            for i = 1:1:4
                for j=2:1:size(imindex,1)
                    if (j == 2)
                        if(imindex(j-1,i)<=12)
                            error1(j,i) = prob(imindex(j-1,i),iter);
                        end
                        for k = 2:1:size(imvalue,1)
                            if(imindex(k,i)>12)
                                error1(j,i) = error1(j,i)+prob(imindex(k,i),iter);
                            end
                        end
                    else
                        if(imindex(j-1,i)>12)
                            error1(j,i) = error1(j-1,i)-prob(imindex(j-1,i),iter);
                        else
                            error1(j,i) = error1(j-1,i)+prob(imindex(j-1,i),iter);
                        end
                    end
                end
                for j = 1:1:size(imvalue,1)
                    if (j<size(imvalue,1))
                        threshold(j+1,i) = (imvalue(j,i) + imvalue(j+1,i))/2;
                    end
                    sign(j,i) = 1;
                    if (error1(j,i)>0.5)
                        error1(j,i) = 1 - error1(j,i);
                        sign(j,i) = -1;
                    end
                end
                threshold(1,i) = imvalue(1,i) - realmin;
            end
            error1(1,:) = 0.5;
            errormin = min(error1);
            final_error = min(errormin);
            [r , c] = find(error1 == final_error);
            test(iter,1) = r(1);
            test(iter,2) = c(1);
            signf = sign(r(1),c(1));
            sign_final(iter) = signf;
            alpha(iter) = 0.5*log((1-final_error)/final_error);
            if (signf == 1)
                for j = 1:1:r(1)-1
                    h(imindex(j,c(1)),iter) = -1;
                end
            else
                for j = r(1):1:size(wtrain,1)
                    h(imindex(j,c(1)),iter) = -1;
                end
            end
            for i=1:1:size(wtrain,1)
                prob(i,iter+1) = prob(i,iter)*exp(-1*alpha(iter)*y(i)*h(i,iter));
            end
            prob(:,iter+1) = (1/sum(prob(:,iter+1))).*prob(:,iter+1);
            iter = iter+1;
        end
        
        % H: Final classifier as a sum up individual classifiers
        H = h*alpha;
        for i = 1:1:length(H)
            if (H(i) <= 0)
                H(i) = -1;
            else
                H(i) = 1;
            end
            if H(i) ~= y(i)
                totalerror = totalerror + 1;
            end
        end
        threshold_value{class} = threshold;
        test_value{class} = test;
        sign_final_value{class} = sign_final;
        class = class + 1;
    end
end
%% Testing
% accuracy computes the name of correctly classified testing trials
accuracy = 0;
% Cell that stores the result of testing for each stimulus
result_value = cell(50,1);
for r = 1:1:50
    test_data = [psth_test1{r}];
    wftest = U'*test_data;
    wftest = wftest';
    result = zeros(50,numtest);
    class = 1;
    count = 0;
    for p = 1:1:49
        for q = p+1:1:50
            threshold = threshold_value{class};
            test = test_value{class};
            sign_final = sign_final_value{class};
            for i = 1:1:size(test_data,2)
                for j = 1:1:maxiter
                    output = 1;
                    if (sign_final(j)  == -1)
                        if (threshold(test(j,1),test(j,2)) < wftest(i,test(j,2)))
                            output = -1;
                        end
                    else
                        if (threshold(test(j,1),test(j,2)) > wftest(i,test(j,2)))
                            output = -1;
                        end
                    end
                    weight_face(i,j) = alpha(j)*output;
                end
                fweight(i) = sum(weight_face(i,:));
                if (fweight(i) > 0)
                    result(p,i) = result(p,i) + 1;
                else
                    result(q,i) = result(q,i) + 1;
                end
                count = count + 1;
            end
            class = class + 1;
        end
    end
    temp = max(result);
    for s = 1:1:length(temp)
        [x1, y1] = find(result == temp(s));
        [x2, y2] = find(x1 == r);
        accuracy = accuracy + length(x2);
    end
    result_value{r} = result;
end
