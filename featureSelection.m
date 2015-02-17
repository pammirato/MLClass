


function [] = featureSelection()
    numSamples = 10000;
    numRelevantFeatures = 3;
    numTotalFeatures = 7;
    mirrorRelevant =0;
    mirrorNonRelevant= 0;

    trainingData = generateData(numSamples,numRelevantFeatures,numTotalFeatures,mirrorRelevant,mirrorNonRelevant);
    
    testData = generateData(numSamples,numRelevantFeatures,numTotalFeatures,mirrorRelevant,mirrorNonRelevant);
    
    %bruteForce(trainingData,testData);
    lassoRegression(trainingData,testData)
    
end




function [errors] = bruteForce(trainingData,testData)

    numFeatures = size(trainingData,2)-1;%-1 because the class is the last column
    featureIndices = [1:numFeatures];
    
    %each cell will be a vector of errors, one for each subset of that size
    errors = cell(1,numFeatures);
    
    %the ith cell will hold all combinations of indices for a subsrt of
    %size i features. i.e. cell 2 = [1 2;1 3;1 4...2 3;2 4...]
    featureIndicesToUse = {};
    
    for i=1:numFeatures
        featureIndicesToUse{i} = combnk(featureIndices,i);
    end
    
    %doLinearRegression(trainingData,testData,featureIndicesToUse{1}(1));
    
    %i goes from 1 to numFeatures
    for i=1:length(featureIndicesToUse)
        subsets =  featureIndicesToUse{i};
        curError= -ones(length(subsets),1);    
        for j=1:length(subsets)
             curError(j) =doLinearRegression(trainingData,testData,subsets(j));
        end
        errors{i} = curError;
    end%for i    
end%end brute force


function [percentWrong] = doLinearRegression(trainingData,testData,featureIndicesToUse)

    %last column in trainingData is the class(y)
    lastCol = size(trainingData,2);
    if(lastCol ~= size(testData,2))
        disp('Error in linear regression data set sizes');
        percentWrong = -1;
        return
    end
    
    
    
    %fitlm does a bias, doesn't need(or want) a column of all 1's
    w = fitlm(trainingData(:,featureIndicesToUse),trainingData(:,lastCol));%.Coefficients{:,1};
    
    w = w.Coefficients{:,1};
    %dont work cause bias
    
    hypothesis = [ones(size(testData,1),1) testData(:,featureIndicesToUse)] *w ;
    
    error = (hypothesis - testData(:,lastCol)).^2;
    
    for i=1:length(hypothesis)
        if(hypothesis(i) < 0)
            hypothesis(i) = -1;
        else
            hypothesis(i) =1;
        end
    end%end for
    
    numWrong = sum(abs(hypothesis-testData(:,lastCol)))/2;
    percentWrong = numWrong/size(testData,1);

end%doLinearRegression


function [errors] = lassoRegression(trainingData,testData)

    numFeatures = size(trainingData,2)-1;%-1 because the class is the last column
    
    %b has coeffieciets for each value of lambda used
    [B STATS] = lasso(trainingData(:,1:end-1),trainingData(:,end),'NumLambda',1000,'LambdaRatio',.0001);
    
    %number of featrues used for each lambda. 
    %numFeatures...numFeatures-1......0
    featuresPerLambda = STATS.DF; 
    lambdas = STATS.Lambda; %all the lambdas used
    
    
    %each cell is a subset of numFeaturesPerLamba
    %cell 1 has all the 1's from numFeaturesLambda, and so on
    
    
    %how many trials we've looked at so far
    soFar = 0;
    %use the largest lambda for each subset size (subset of features)
    trialsToUse = zeros(numFeatures,1);
    for i=1:numFeatures
               %gets number of lambdas used per subset size
        trialsDone = length(featuresPerLambda(featuresPerLambda == i));
        if(trialsDone >0)
            trialsToUse(i) =  trialsDone+ soFar;
            soFar = trialsDone + soFar;
        else
            trialsToUse(i) = -1;
        end
    end
    %so for each possible number of features to use,  pick one
    %set of coefficients coresponding to the largest lambda for that number
    
    
   %get errors
   errors = -ones(numFeatures,1);
   for i=1:numFeatures
       if(trialsToUse(i) <1)
           continue;
       end;
       errors(i) = getError(B(1:numFeatures,trialsToUse(i)),testData);
   end
end%end doLassoRegression





function [data] =  generateData(numSamples,numRelevantFeatures,numTotalFeatures,mirrorRelevant,mirrorNonRelevant)

    %generate the relevant features for class 1
               %mean -4, standard deviation 2
    xPos =  -4 + 2.*randn(numSamples,numRelevantFeatures) + (2 + 4.*randn(numSamples,numRelevantFeatures));
    xNeg =  4 + 2*randn(numSamples,numRelevantFeatures) +  (-2 + 4*randn(numSamples,numRelevantFeatures));
    
    %generate the non relevant features, just random numbers
    nonRel = -6 + (6--6).*rand(numSamples*2,numTotalFeatures-numRelevantFeatures);
    
    %put half in postive class, half in negative class
    xPos = [xPos nonRel(1:numSamples,:)];
    xNeg = [xNeg nonRel(numSamples+1:size(nonRel,1),:)];
    
    if(mirrorRelevant)
        if(numTotalFeatures > numRelevantFeatures)
            %make last feature(non-relevant) = -first feature(relevant)
            xPos(:,size(xPos,2)) = -xPos(:,1);
        else
            disp('no non relevant feature to make a mirror');
        end
    end%if mirror relevant

    if(mirrorRelevant)
        if(numTotalFeatures > (numRelevantFeatures+1))
            %make last feature = -next to last feature(both are not relevant
            xPos(:,size(xPos,2)) = -xPos(:,size(xPos,2)-1);
        else
            disp('not enough non relevant feature to make a mirror');
        end
    end%if mirror relevant
    
    %add classification
    xPos = [xPos ones(size(xPos,1),1)];
    xNeg = [xNeg -ones(size(xNeg,1),1)];
    
    data = [xPos;xNeg];
end




function [percentError] = getError(hypothesisCoeff,testData)

    hypothesis = testData(:,1:end-1) *hypothesisCoeff ;
    
    for i=1:length(hypothesis)
        if(hypothesis(i) < 0)
            hypothesis(i) = -1;
        else
            hypothesis(i) =1;
        end
    end%end for
    
    numWrong = sum(abs(hypothesis-testData(:,end)))/2;
    percentError = numWrong/size(testData,1);
end