function [startPosValid,lengthsValid,thresholdPos,thresholdNeg] = getPosition(data,probability)

    % Calculate statistical threshold for the current signal
    dataDiff = diff(data); 
    
    posChance = length(find(dataDiff>0))/length(dataDiff);
    negChance = 1-posChance;
        
    thresholdPos = round(-probability/log2(posChance));
    thresholdNeg = round(-probability/log2(negChance));
    
    
    signVec = sign(dataDiff);
    signVec_corrected = signVec;
    % In some experiments we noticed we were missing out on some obvious
    % activity segments because a single point in the segment was
    % "misalinged" due to measurement noise. To counter this while still
    % maintaing some measure of stringency, we decided that in cases where
    % a center point in a stretch of 5 points is misalinged with the rest
    % of the points (e.g., a decrease in an otherwise increasing environment),
    % we will manually correct this.
    for i = 3:length(dataDiff)-2
        if (signVec(i-2) == signVec(i-1)) && (signVec(i-1) == signVec(i+1)) ...
                && (signVec(i+1) == signVec(i+2)) && (signVec(i-1)~=signVec(i))
            signVec_corrected(i) = signVec(i-1);
        end
    end
    
    signVec_corrected(signVec_corrected==-1) = 0;
    % runLengthEncode function looks for contigous segments of the same
    % value in a vector
    [startPos,lengths, ~] = runLengthEncode(signVec_corrected');
    points = lengths > min([thresholdPos,thresholdNeg]); % Initial filterring of points
    lengthsValid = lengths(points);
    startPosValid = startPos(points);
end

