function [Kinder,MonitoringMatrix] = swap_mutation(Kinder,MutationWhk,p,n)


MonitoringMatrix = cell(p,1);


for z1=1:p
    index = 1;
    for z2=1:n
        kappa = rand(1);
        if kappa < MutationWhk
            tausch_gen = round(n*rand(1));
            while tausch_gen == z2 || tausch_gen == 0
                tausch_gen = round(n*rand(1));
            end
            tmp_gen = Kinder{z1,1}(1,z2);                        % Gen 1 zwischenspeichern
            Kinder{z1,1}(1,z2) = Kinder{z1,1}(1,tausch_gen);          % Gen 2 in Gen 1 speichern
            Kinder{z1,1}(1,tausch_gen) = tmp_gen;                % Gen 1 (tmp) in Gen 2 speichern
            MonitoringMatrix{z1,1}(index,1) = z2;
            MonitoringMatrix{z1,1}(index,2) = tausch_gen;
            index = index + 1;
        end
    end
end