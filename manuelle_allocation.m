function [Population] = manuelle_allocation(Population,p,geneM,geneR,permutation)

% Allocation
    % Mensch
if geneM(1,1) ~= 0
    for z1 = 1:p
        for s = 1:length(geneM)
            Population{z1,1}(2,geneM(1,s)) = 1;
            Population{z1,1}(3,geneM(1,s)) = 0;
        end
    end
end

    % Roboter
if geneR(1,1) ~= 0
    for z1 = 1:p
        for s = 1:length(geneR)
            Population{z1,1}(2,geneR(1,s)) = 0;
            Population{z1,1}(3,geneR(1,s)) = 1;
        end
    end 
end      


% Permutation
if permutation(1,1) ~= 0
    [~,spalten] = size(permutation);
    for z1=1:p
        for s = 1:spalten
            [~,x] = find(Population{z1,1}(1,:) == permutation(2,s));        % Operation x finden, die die neue Permu in Zeile 2 bisher besitzt
            Population{z1,1}(1,x) = Population{z1,1}(1,permutation(1,s));   % Operation x erhält alte Permutation 
            Population{z1,1}(1,permutation(1,s)) = permutation(2,s);        % Neue Permutation an Operation geben
        end
    end
end