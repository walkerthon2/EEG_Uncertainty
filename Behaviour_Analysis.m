%Function to extract behavioural data
%Start defines the first participant to analyse, finish the final
%participant. 
function [output] = Behaviour_Analysis(start, finish)

close all
format bank
clear rtData
clear ctData

%Basic exp info
numTrials = 512;
numBlocks = 32;

%Divided by two for certain/uncertain
trialsPerBlock = numTrials / numBlocks / 2;

%Defining a struct for certain and uncertain data. 
rtData.cert = nan(finish-start+1, numBlocks/2);
ctData.cert = nan(finish-start+1, numBlocks/2);

rtData.unc = nan(finish-start+1, numBlocks/2);
ctData.unc = nan(finish-start+1, numBlocks/2);

%Defining a struct for output
output.rtCert = NaN(finish-start+1, numBlocks/4);
output.rtUnc = NaN(finish-start+1, numBlocks/4);
output.ctCert = NaN(finish-start+1, numBlocks/4);
output.rtUnc = NaN(finish-start+1, numBlocks/4);

for s = start:finish
    
    %Load Data in for each Subject
    subCondName = ['C1_S' int2str(s), '_ET0'];
    dataLoc = ['Data\', subCondName];
    load(dataLoc, 'DATA');
    
    %Split data into cert and uncert
    certData = DATA.part1(DATA.part1(:, 6) == 1, :);
    uncData = DATA.part1(DATA.part1(:, 6) == 2, :);
    
    %Calculate Reaction Times for certain and uncertain
    rtData.cert(s, 1:numBlocks) = dataCondense(certData, 13, numBlocks, 13, trialsPerBlock);
    rtData.unc(s, 1:numBlocks) = dataCondense(uncData, 13, numBlocks, 13, trialsPerBlock);

    %Calculate Performance.
    ctData.cert(s, 1:numBlocks) = dataCondense(certData, 8, numBlocks, 13, trialsPerBlock);
    ctData.unc(s, 1:numBlocks) = dataCondense(uncData, 8, numBlocks, 13, trialsPerBlock);
    
    %Collapse for output into epochs of 32 trials
    output.rtCert(s, 1:numBlocks / 4) = blockCondense(4, numBlocks, rtData.cert(s, :));
    output.rtUnc(s, 1:numBlocks / 4) = blockCondense(4, numBlocks, rtData.unc(s, :));
    output.ctCert(s, 1:numBlocks / 4) = blockCondense(4, numBlocks, ctData.cert(s, :));
    output.ctUnc(s, 1:numBlocks / 4) = blockCondense(4, numBlocks, ctData.unc(s, :));
    
    
end

Plot(1, 1:8) = nanmean(output.ctCert, 1);
Plot(2, 1:8) = nanmean(output.ctUnc, 1);

plot(Plot'); ylim([.4 1]); legend('Certain', 'Uncertain');

%Exclusions if you want them
output.ctCert(2, :) = NaN;
output.ctCert(3, :) = NaN;
output.ctCert(5, :) = NaN;
output.ctCert(12, :) = NaN;

output.ctUnc(2, :) = NaN;
output.ctUnc(3, :) = NaN;
output.ctUnc(5, :) = NaN;
output.ctUnc(12, :) = NaN;

figure(2)
Plot(1, 1:8) = nanmean(output.ctCert, 1);
Plot(2, 1:8) = nanmean(output.ctUnc, 1);

plot(Plot'); ylim([.4 1]); legend('Certain', 'Uncertain');

end

function dataFix = dataCondense(Data, column, blockTotal, reactionTime, trialsPerBlock)

dataFix = nan(1,blockTotal);

%trialsPerBlock = 32;

for b = 1:blockTotal                                                       %Read in each block
    %Create logical array for each block for the relevant block    
    stdTime = std(Data((trialsPerBlock * (b - 1) + 1):trialsPerBlock * b, reactionTime));
    meanTime = mean(Data((trialsPerBlock * (b - 1) + 1):trialsPerBlock * b, reactionTime));
    
    %Removes trials with times over or under 2 std deviations from the mean
    %on that block
    test = Data(:, 2) == b & Data(:, reactionTime) < meanTime + stdTime * 2 ...
    & Data(:, reactionTime) > meanTime - stdTime * 2; %Select block and remove outlier reaction times                           
    %end
    
    if ismember(1, test) == 1                                              %Test to see if that pattern was present in the block
        dataTemp = Data(test, column);                                      %Write a temporary array for all of the participant's reaction time's in a block
        test = double(test);    
        for i = 1:trialsPerBlock
            if Data((b - 1) * trialsPerBlock + i, reactionTime) > meanTime + stdTime * 2 ...
            || Data((b - 1) * trialsPerBlock + i, reactionTime) < meanTime - stdTime * 2
                test(i, 1) = nan;
            end
        end
        dataFix(:, b) = nanmean(dataTemp);                                  %Condense those results down to 1 mean value
        %dataFix(:, b) = dataTemp;   
    else
        dataFix(:, b) = NaN;                                                %If there's no values in that block, return a NaN
    end       
        
end
    
end


function condensedData = blockCondense(block, blockTotal, data)             %Call in the amount to condense, total amount of blocks, and data to condense
blockCounterStart = 1;                                                      %Set counters for position in data
blockCounterEnd = block;
for c = 1:blockTotal/block
             
    tempBlock = data(:, blockCounterStart:blockCounterEnd);

    tempData(:, c) = nanmean(tempBlock, 2);
        
    blockCounterStart = blockCounterStart + block;
    blockCounterEnd = blockCounterEnd + block;
end
    
    
condensedData = tempData;
end