
%main function
function [] = hw1()

    %generate points with label +1 and -1
    posPoints = generatePointsBeta(10,100,0,40,1);
    negPoints = generatePointsBeta(30,100,0,40,-1);
    
    %put all the points in a single matrix
    allPoints = horzcat(posPoints,negPoints)

    
    %plot the 2D points
    plotPoints(posPoints,'b.')
    hold on
    plotPoints(negPoints,'r.')
    hold on
    
    
    %plot the loss functions
    plotLoss(allPoints)
    
    
    
end

%plots the 2d points

%points array has the format     x1, x2, ...
%                                y1, y2, ...
%                                label1, label2, ...
function [] = plotPoints(points,format)

    plot(points(1,:),points(2,:),format)

end


%gnerate points with each points x,y coordinates chosen from
%a beta distribution with the given mean.
function [points] = generatePointsBeta(mean,number,lowRange,highRange, sign)

    % make the mean in [0,1] for the beta dist.
    betaMean = mean / (highRange-lowRange) ;

    %equation for mean of beta dist rearranged - Wikipedia
    %to get the alpha,beta parameters
    baRatio = (1/betaMean) - 1;
    
    
    beta = 1; %just because  :(
    alpha = beta/baRatio;
    
    %generate the random points
    points = random('beta',alpha,beta,[2,number]);
    
    %scale the points from [0,1] to the given range
    points = (points*highRange) + lowRange;
    
    signs(1:number) = sign;
    
    points = [points;signs];

end



%plots the zero one and hinge loss for the points using the function
%giveen in the homework description
function [] = plotLoss(points)


    %range of bias to check
    b =0;
    
    %holds the value of loss for each bias checked
    zolb(1:length(b)) = 0;
    hilb(1:length(b)) = 0;

    %let w1,w2 range from -10 to 10 by 1.
    w1 = -10:.5:10;
    w2 = -10:.5:10;
    
    %hold the loss value
    zoLoss  = zeros(length(w1),length(w2));
    hLoss  = zeros(length(w1),length(w2));
    
    %zoLoss  = zeros(360);
    %hLoss  = zeros(360);
    
    
    numPoints = length(w1);
    
    
    
    for i=1:numPoints,
        for j=1:numPoints
     
     %for angle=0:1:360 , i didnt really get the angle thing
            
            for k=1:length(zolb),
                zolb(k) = zeroOneLoss(w1(i),w2(j),b(k),points);
                hilb(k) = hingeLoss(w1(i),w2(j),b(k),points);
                %zolb(k) = zeroOneLoss(cos(angle),sin(angle),b(k),points);
                %hilb(k) = hingeLoss(cos(angle),sin(angle),b(k),points);
            end
            %compute the zero one loss for this choice of w1,w2
            
            
            zoLoss(i,j) = min(zolb);
            hLoss(i,j) = min(hilb);
            
            
            
    % end
        end
        
    end
    
    
    figure
    [X,Y] = meshgrid(w1,w2);
    surf(X,Y,zoLoss)
    
    figure
    surf(X,Y,hLoss)

end




function [loss] = zeroOneLoss(w1,w2,b, points)

    idontknowmatlab = size(points);
    
    loss = 0;
    
    
    for i=1:idontknowmatlab(2),
        
        
        %value of the decision function
        f = sign(w1*points(1,i) + w2*points(2,i)+b);
        
        %actual label
        y = points(3,i);
        
        loss = loss + .5*(1-f*y);
    end
end



function [loss] = hingeLoss(w1,w2,b, points)

    idontknowmatlab = size(points);
    
    loss = 0;
    
    
    for i=1:idontknowmatlab(2),
        %              max(0,1-      f(x)                                  *y)
        loss = loss + max(0, 1 - (w1*points(1,i) + w2*points(2,i) + b)*points(3,i));   
    end
end
