function [Kinder,MonitoringMatrix] = ein_bit_mutation (Kinder,MutationWhk,p,n)

MonitoringMatrix = cell(p,1);

for z1=1:p
    index = 1;
    for z2=1:n
        kappa = rand(1);
        if kappa < MutationWhk
            if Kinder{z1,1}(2,z2) == 1
               Kinder{z1,1}(2,z2) = 0;
            else
                Kinder{z1,1}(2,z2) = 1;
            end
            MonitoringMatrix{z1,1}(index,1) = z2;
            index = index + 1;
        end
    end
end

for z1 = 1:p
    Kinder{z1,1}(3,:) = ones(1,n) - Kinder{z1,1}(2,:);
end
