function Kinder = mixed_mutation(Kinder,p,n,genom,mut_whk)
% Der Parameter Genom ermöglicht die Anwendung auf der mixed Mutation auf Sequenz (Genom=1) und Zuordnung(Genom=2)

for z1 = 1:p
    random = rand(1);
    
    if random < mut_whk
    
        schnitt1 = round(n*rand(1));
        while schnitt1 == 0
            schnitt1 = round(n*rand(1));
        end

        schnitt2 = round(n*rand(1));
        while schnitt1 == schnitt2 || schnitt2 == 0
            schnitt2 = round(n*rand(1));
        end

        if schnitt2 < schnitt1
            [schnitt1,schnitt2] = deal(schnitt2,schnitt1);
        end

        abschnitt = Kinder{z1,1}(genom,schnitt1:schnitt2);

        [~, nRow] = size(abschnitt);
        abschnitt(1, :) = abschnitt(1, randperm(nRow));

        Kinder{z1,1}(genom,schnitt1:schnitt2) = abschnitt;
    end
    
end