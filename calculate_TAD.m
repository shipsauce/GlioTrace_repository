function TAD = calculate_TAD(xs, ys)
% This function returns the turning angle distribution of a sequence by 
% calculating the relative turning angles in a sequence of consecutive coordinates.
%
% Input parameters:
% xs - x-coordinates
% ys - y-coordinates
%
% Output parameters:
% TAD - vector of turning angles
%
% @authors: Madeleine Skeppås
% @date: 16012025

theta = [];

% Calculate x,y-displacement in consecutive frames
v = [xs(2:end) - xs(1:end-1), ys(2:end) - ys(1:end-1)];

if(~isempty(v))
    v2 = v(or((~v(:,1)==0), (~v(:,2)==0)),:); % Ignore steps without displacement

    for j=1:height(v2)-1
        % Calculate the cosine of the relative turning angle
        exp = dot(v2(j,:),v2(j+1,:)) / (norm(v2(j,:)) * norm(v2(j+1,:)));
        
        % Convert to radians and handle edge cases
        if(exp <= -1)
            theta(j) = acos(-1); % pi radians/ 180 degrees turn
        elseif(exp >= 1)
            theta(j) = acos(1); % 0 degrees
        else
            theta(j) = acos(dot(v2(j,:),v2(j+1,:)) / (norm(v2(j,:)) * norm(v2(j+1,:))));
        end
    end
end

TAD = theta(find(theta));

end