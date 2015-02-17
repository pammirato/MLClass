


function [] = hw3()

    numValSets = 5;
    numTestingPoints = 100;
    
    kNNRange = 10;
    
    
    %numTrials = 10;
    
    numTrainingPoints = 10:10:100;
    numValidationPoints = numTrainingPoints/numValSets;
    avgValErrorsKNN = zeros(size(numTrainingPoints,2),kNNRange);
    testErrorsKNN = zeros(size(numTrainingPoints,2),kNNRange);
    trainingErrorsKNN = zeros(size(numTrainingPoints,2),kNNRange);
    
    avgValErrorsLR = zeros(size(numTrainingPoints,2));
    testErrorsLR = zeros(size(numTrainingPoints,2));
    trainingErrorsLR = zeros(size(numTrainingPoints,2));

    %create two functions, one foreach class, -1 or 1
    %each is a mixture or Gaussians
    f1 = inline('normpdf(x,1,2) + normpdf(x,-6,1)','x'); %class 1
    f2 = inline('normpdf(x,6,1) + normpdf(x,-1,2)','x'); %class -1
    
    
    
    for trial=1:size(numTrainingPoints,2)
         %init training sets, each row has x,y,class
        f1Train = ones(3,numTrainingPoints(trial));
        f2Train(1:3,1:numTrainingPoints(trial)) = -1;


        %generate the points
        f1Train(1,:) = generateValues(f1, size(f1Train,2));
        f1Train(2,:) = generateValues(f1, size(f1Train,2));
        f2Train(1,:) = generateValues(f2, size(f2Train,2));
        f2Train(2,:) = generateValues(f2, size(f2Train,2));






        %init testing sets, each row has x,y,class
        f1Test = ones(3,numTestingPoints);
        f2Test(1:3,1:numTestingPoints) = -1;

        %generate testingPoints
        f1Test(1,:) = generateValues(f1, size(f1Test,2));
        f1Test(2,:) = generateValues(f1, size(f1Test,2));
        f2Test(1,:) = generateValues(f2, size(f2Test,2));
        f2Test(2,:) = generateValues(f2, size(f2Test,2));

        testingSet = [f1Test f2Test];


        plot(f1Train(1,:),f1Train(2,:),'.b');
        hold on;
        plot(f2Train(1,:),f2Train(2,:),'.r');
        hold on;
        title('Last Set of Training Points');


        for k=1:10

            totalValError = 0;
            for i=1:numValSets
                valStart = (i-1)*numValidationPoints(trial) +1;
                valEnd = valStart+numValidationPoints(trial)-1;

                valSet = [f1Train(1:3,valStart:valEnd) f2Train(1:3,valStart:valEnd)];

                trainingSet = setdiff( [f1Train f2Train]', valSet', 'rows')'; 

                totalValError = totalValError + errorKNN(trainingSet,valSet,k);

            end%for i

            avgValErrorsKNN(trial,k) = totalValError/numValSets;
            
            %now with full trainingSet
            fullTrainingSet = [f1Train f2Train];

            trainingErrorsKNN(trial,k) = errorKNN(fullTrainingSet,fullTrainingSet,k);

            testErrorsKNN(trial,k) = errorKNN(fullTrainingSet, testingSet, k);

            %fprintf('K = %d ,  ValE = %d  ,  TestE = %d   ,  TrainE = %d \n',k,averageValError,testError, trainingErros(trial);

        end%for k=1:10
        
        %some linear regression for fun
        
        totalValError = 0;
        for i=1:numValSets
            valStart = (i-1)*numValidationPoints(trial) +1;
            valEnd = valStart+numValidationPoints(trial)-1;

            valSet = [f1Train(1:3,valStart:valEnd) f2Train(1:3,valStart:valEnd)];

            trainingSet = setdiff( [f1Train f2Train]', valSet', 'rows')';
            w = linearRegression(trainingSet);

            totalValError = totalValError + errorLR(w,valSet);

        end%for i

        avgValErrorsLR(trial,k) = totalValError/numValSets;

        %now with the full training set
        fullTrainingSet = [f1Train f2Train];
        w = linearRegression(fullTrainingSet);
        trainingErrorsLR(trial) = errorLR(w,fullTrainingSet);

        testErrorsLR(trial) = errorLR(w,testingSet);
    
    end%for trials
    
    
   
    
    
    
    for k=1:10
        
        figure;
        plot(numTrainingPoints,trainingErrorsKNN(:,k),'.--b', ...
            numTrainingPoints,testErrorsKNN(:,k),'x-r', ...
            numTrainingPoints,avgValErrorsKNN(:,k),'+:g');
        legend('TrainingError','Test Error', 'Cross Val Avg Error');
        xlabel('Training Set Size(including 1/5 validation)');
        ylabel('Percent Error');
        titlee = sprintf('K = %d',k);
        title(titlee);
        
    end%for k=1:10   the second one
    
    
    figure;
    plot(numTrainingPoints,trainingErrorsLR,'.--b', ...
        numTrainingPoints,testErrorsLR,'x-r', ...
        numTrainingPoints,avgValErrorsLR,'+:g');
    legend('TrainingError','Test Error', 'Cross Val Avg Error');
    xlabel('Training Set Size(including 1/5 validation)');
    ylabel('Percent Error');
    titlee = 'Linear Regression';
    title(titlee);
    

end %hw3








%generates n values from the pdf f
function [samples] = generateValues(f,n)

    samples = zeros(1,n);

    %seed the uniform random number generator
    rand('seed',now);
    %a vector that has the range of the experiment
    x = -10:.1:10;
    
  
    
    %create Gaussian for rejection sampling
    q = inline('normpdf(x,0,4)','x');
    %q2 = inline('normpdf(x,0,4)','x');
    
    %get ratios of scalingg for so q's contain f's
    c = max(f(x)./q(x));
    %c2 = max(f2(x)./q2(x));
    
    
    count = 0;
    
    while (count < n)
        %draw potential sample from q1,q2
        s = normrnd(0,4);
        %s2 = normrnd(0,4);


        %get likelihood this point would be in f's
        ls = f(s) / (c*q(s));
        %ls2 = f2(s2) / (c*g2(s2));

        %generate a random number uniformally in (0,1)
        u = rand();

        if(ls > u)
            count = count +1;
            samples(count) = s;
        end%uf ls>u
        
    end %while count
end%generate Points




function [w] = linearRegression(trainingSet)
    
    %add the 1 for w0
    trainingSet = [ones(1,size(trainingSet,2)) ; trainingSet];


    w = pinv(trainingSet(1:3,:)') * trainingSet(4,:)';    
end


function [percentWrong] = errorLR(w, testingSet)

    %add for w0
    testingSet = [ones(1,size(testingSet,2)) ; testingSet];

    totalWrong = 0;
    for i=1:size(testingSet,2)
        class = w * testingSet(1:3,i)';
        if(class < 0)
            class = -1;
        else
            class = 1;
        end    
        if(class ~= testingSet(4,i))
            totalWrong = totalWrong +1;
        end
    end
    
    percentWrong = totalWrong/size(testingSet,2);
end%errorLR


function [class] = findKNNDumb(point, pointSet, k)

    count = 0;
    
    %holds the k nearest points and their distances and classes
    knns = zeros(4,k);

    
    %distance to fisrt point
    maxNearestDist = ( point(1,1) - pointSet(1,1) )^2 + ( point(2,1) - pointSet(2,1) )^2;
    
    
    for i=1:size(pointSet,2)
        
        dist = (point(1,1) - pointSet(1,i))^2 + (point(2,1) - pointSet(2,i))^2;
        
        %if this point is as least as close as the farthest nn
        if(dist <= maxNearestDist || count < k)
            %if we havnet found knn yet, just add it
            if(count < k)
                count = count + 1;
                knns(1:4,count) = [dist ; pointSet(1:3,i)];
                if(dist > maxNearestDist)
                    maxNearestDist = dist;
                end
            
            else
                knns = sortrows(knns');%sort the knns by distance
                knns = knns';
                if(dist < knns(1,k))
                    knns(1:4,k) = [dist ; pointSet(1:3,i)]; %replace the farthest point 
                    if(k >1 && dist < knns(1,k-1))
                        maxNearestDist = knns(1,k-1);
                    else
                        maxNearestDist = dist;
                    end
                end
                
            end
        
        end
        
    end %for
    
    %sum up the classification of all neighbors
    class = sum(knns(4,:));
    
    if(class < 0)
        class = -1;
    else
        class = 1;
    end
    

end%findknnndumb


%returns percentage of points not classified correctly
function [percentWrong] = errorKNN(trainingSet, testingSet, k)

    totalWrong = 0;

    for i=1:size(testingSet,2)
        
        class = findKNNDumb( testingSet(1:2,i), trainingSet, k);
        
        if(class ~= testingSet(3,i))
            totalWrong = totalWrong +1;
        end
        
    end%for i

    percentWrong = totalWrong/size(testingSet,2);

end

