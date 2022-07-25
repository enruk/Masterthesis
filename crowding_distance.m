function [Population] = crowding_distance (Population,p)

    % Genotypisch
    % Hamming - Abstand berechnen
for z1 = 1:p
    Population{z1,7} = zeros(p,2);
end


for z1 = 1:p
    for z2 = 1:p
        HammingMatrix = [Population{z1,1}(2,:);Population{z2,1}(2,:)];
        HammingDistanz = size(HammingMatrix,2)*pdist(HammingMatrix,'Hamming');
        Population{z1,7}(z2,1) = HammingDistanz;
    end
end


    % Phänotypisch
    % durchschnittliche Zeitdistanz berechnen

for z1 = 1:p
    for z2 = 1:p
        Population{z1,7}(z2,2) = mean(abs(Population{z1,1}(4,:) - Population{z2,1}(4,:)));
    end
end