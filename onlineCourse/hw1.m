

function [] = hw1()

    display = 1;
    totalIterations = 0;
    numTrials = 1;%000;
    numTrainingPoints = 10;
    
    for trials=1:numTrials

        target = getRandomLine(); 
        %clf;
        %plot(target(1,:),target(2,:),'r')

        trainingSet = rand(3,numTrainingPoints)*2 -1;

        %classify training data
        for i=1:length(trainingSet)
            trainingSet(3,i) = targetValue(target, trainingSet(1:2,i));
        end


        %save the target values, these will get overwritten in the trainingSet
        targetValues = trainingSet(3,:);


        %add the 1 for w0
        trainingSet = [ones(1,size(trainingSet,2)) ; trainingSet];



        %labels = cellstr(num2str(targetValues'));
        %hold on
        %plot(trainingSet(2,:), trainingSet(3,:),'.b')
        %text(trainingSet(2,:), trainingSet(3,:), labels, 'VerticalAlignment','bottom', ...
        %                         'HorizontalAlignment','right')

        if(display)
            plotTargetAndPoints(target,trainingSet,targetValues);
        end
        
        w = [0;0;0];


        done = 0;
        index = -1;%for the chosen point to move the hypothesis around

        numIterations = 0;
        

        while not(done)
            trainingSet = classifyPoints(w,trainingSet);
            numIterations = numIterations +1;

            %plot everything
            if(display)
                plotTargetAndPoints(target,trainingSet,targetValues)
            
            
                %now plot hypothesis
                hold on
                if (w(3) == 0)
                   slope = 0;
                   intercept = 0;
                else
                    slope = (-w(2))/w(3);
                    intercept = -w(1)/w(3);
                end

                xx = -1:1;
                yy = slope*xx + intercept;
                plot(xx,yy,'--g');

                if index > 0
                    plot(trainingSet(2,index), trainingSet(3,index), '.r');
                end

                axis([-1.2,1.2,-1.2,1.2]);
                pause()
            end


            if(sum(trainingSet(end,:) - targetValues) ~= 0)
                index = chooseMisclassifiedPoint(trainingSet,targetValues);



                if(index == -1)
                    disp('Error')
                end

                w = w + targetValues(index)*trainingSet(1:3,index);
            else
                done =1;
                totalIterations = totalIterations + numIterations;
            end


        end %while not done
    end
    
    avgIterations = totalIterations/numTrials
       
    
end% hw1 function




function [] = plotTargetAndPoints(target,trainingSet,targetValues)
    clf;
    tslope = (target(2,2) - target(2,1)) / (target(1,2) - target(1,1));
    tintercept = target(2,1) - tslope*target(1,1);
    xx = -1:1;
    tyy = xx * tslope + tintercept;
    plot(xx,tyy,'r');
    %plot(target(1,:),target(2,:),'r')%plot target function line segment
    hold on
    plot(trainingSet(2,:), trainingSet(3,:),'.b')%plot training set points
    %label the trainingSet points with there target values
    text(trainingSet(2,:), trainingSet(3,:),cellstr(num2str(targetValues')) , 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')
    
    %label with hypothesis, on the right                     
    text(trainingSet(2,:), trainingSet(3,:),cellstr(num2str(trainingSet(end,:)')) , 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','left')                      
                     


end





function [index] = chooseMisclassifiedPoint(points, targetValues)
    index = -1;
    for i=1:size(points,2)
        
        if points(end,i) ~= targetValues(i)
            index = i;
            break
        end
    end
end


function [done] = classifyPoints(w,x)

    %add the 1 for w0
    if(size(x,1) < 3)
        x = [ones(size(x,2)) ; x(1,:) ; x(2,:)] ;
    end
    
    done = hypothesis(w,x);

end


function [class] = hypothesis(w,x)

    if sum(w) == 0
        x(size(x,1),:) = 0; %put 0s in the last row
        class = x;
        return
    end
    
    x(size(x,1),:) = sign(w'*x(1:3,:));
    class = x;
        
    

end

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







%NO
function [x] = randomBinaryVector(x,bLow,bHigh)

    for i=1:length(x),
        if x(i) < .5;
            x(i) = -1;
        else
            x(i) = 1;

        end
    end
end