function y = CalculateAngleDifference(x1,x2,angleType)
% Calkculate the angle difference in degrees

d = abs(x1-x2);
switch angleType
    case 'azi'
        %y = abs(min(360-d,d));
        dc = 360-d;
    case 'ele'
        %y = abs(min(180-d,d));
        dc = 180-d;
    otherwise
        error("Must be 'azi' or 'ele'.")
end

[dmin,idx] = min([d,dc]);
% taking x1 - x2
if idx == 1
    if x1 > x2
        y = dmin;
    else
        y = -dmin;
    end
else
    if x1 > x2
        y = -dmin;
    else
        y = dmin;
    end
end