

function []= hw2()
   
    linearRegression();

end






function [] = linearRegression()


    display = 0;
    totalIterations = 0;
    numTrials = 1000;
    numTrainingPoints = 100;
    numTestingPoints = 1000;
    totalTrainWrong= 0;
    totalTestWrong = 0;
    
    
    for trial=1:numTrials

        target = getRandomLine(); 


        %generate TrainingPoints
        trainingSet = rand(3,numTrainingPoints)*2 -1;

        %classify training data
        for i=1:length(trainingSet)
            trainingSet(3,i) = targetValue(target, trainingSet(1:2,i));
        end


        %save the target values, these will get overwritten in the trainingSet
        targetValues = trainingSet(3,:);
        
        %add the 1 for w0
        trainingSet = [ones(1,size(trainingSet,2)) ; trainingSet];


        w = pinv(trainingSet(1:3,:)') * targetValues';

        
        %how good was training?
        numWrong = 0;

        hypothesisValues = w'*trainingSet(1:3,:);

        for i=1:length(hypothesisValues)

            if(sign(hypothesisValues(i)) ~= targetValues(i))
                numWrong = numWrong +1;
            end
        end%for i

        totalTrainWrong = totalTrainWrong +numWrong;
        
        
        
        
        
        %now testttttttttttttttttttttttttttttttttttt
        
        
        %generate testing set
        testingSet = rand(3,numTestingPoints)*2 -1;
        
        %add the 1 for w0
        testingSet = [ones(1,size(testingSet,2)) ; testingSet];
        
        %classify testing data for comparision
        for i=1:length(testingSet)
            testingSet(4,i) = targetValue(target, testingSet(2:3,i));
        end
        
        
        
        
        
        
        numWrong = 0;

        hypothesisValues = w'*testingSet(1:3,:);

        for i=1:length(hypothesisValues)

            if(sign(hypothesisValues(i)) ~= testingSet(4,i))
                sign(hypothesisValues(i));
                testingSet(4,i);
                numWrong = numWrong +1;
            end
        end%for i

        totalTestWrong = totalTestWrong +numWrong;
        
        

    end%for trial
    
    averageTrainPercentageWrong = totalTrainWrong/numTrials/numTrainingPoints
    averageTestPercentageWrong = totalTestWrong/numTrials/numTestingPoints
    
end



%from hw1

%returns the target value for the point
function [target] = targetValue(target, point)
    %this indexing is ok, the ones have not been added to points
    target = (not( left(target(:,1), target(:,2), point))) *2 -1 ; %gives -1 or +1
end


function [line] = getRandomLine()

    line = [getRandomPoint() ; getRandomPoint()];
    
    %put higher point on top
    if line(2,1) > line(2,2)
        
        temp = line(:,1);
        line(:,1) = line(:,2);
        line(:,2) = temp;
    end
end


function [point] = getRandomPoint()

    point = rand(1,2)*2 -1;
    %point = randomBinaryVector(point,-1,1);
end



%returns 1 if c is left of directed line from a to b. 
function [tf] = left(a,b,c)

    tf = area2(a,b,c) > 0;
end


function [area] = area2(a,b,c)
    
    area = (b(1) - a(1)) * (c(2)-a(2)) - (c(1)-a(1)) * (b(2)-a(2));
end


%from hw1




















%for number 1 and 2
function [averages] = coinFlipExperiment()

    numCoins = 1000;
    numFlips = 10;

    numExperiments = 100000;
    frlFractioins = zeros(numExperiments,3);

    for i=1:numExperiments

        results = floor(rand(numFlips,numCoins)*2);

        firstCoin = results(:,1);
        randomCoin = results(:,floor(rand()*numCoins)+1);
        lowHeadsCoinValue = min(sum(results,1));

        fractionFirst = sum(firstCoin)/numFlips;
        fractionRand = sum(randomCoin)/numFlips;
        fractionLow = lowHeadsCoinValue/numFlips;

        frlFractions(i,:) = [ fractionFirst, fractionRand, fractionLow];

    end


    averages = sum(frlFractions,1)/numExperiments

end
