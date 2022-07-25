function [Kinder] = two_point_crossover(Kinder,Paare,p,n,AlloCrossWhk)

for z1 = 1:p
    
    random = rand(1);
    if random < AlloCrossWhk
        
        Kinder{z1,1}(2,:) = zeros(1,n);

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


        if schnitt1 == 1 && schnitt2 < n        % Fall 1: Abschnitt beginnt bei 0
            Kinder{z1,1}(2,schnitt1:schnitt2)= Paare{z1,1}(2,schnitt1:schnitt2);
            Kinder{z1,1}(2,schnitt2+1:n)= Paare{z1,2}(2,schnitt2+1:n);


        elseif schnitt1 == 1 && schnitt2 == n   % Fall 2: Abschnitt ist kompletter String
            Kinder{z1,1}(2,schnitt1:schnitt2)= Paare{z1,1}(2,schnitt1:schnitt2);


        elseif schnitt1 > 1 && schnitt2 == n    % Fall 3: Abschnitt endet mit n
            Kinder{z1,1}(2,schnitt1:schnitt2) = Paare{z1,1}(2,schnitt1:schnitt2);
            Kinder{z1,1}(2,1:schnitt1-1) = Paare{z1,2}(2,1:schnitt1-1);


        else                                    % Fall 4: Abschnitt irgendwo in der Mitte
            Kinder{z1,1}(2,schnitt1:schnitt2) = Paare{z1,1}(2,schnitt1:schnitt2);
            Kinder{z1,1}(2,1:schnitt1-1) = Paare{z1,2}(2,1:schnitt1-1);
            Kinder{z1,1}(2,schnitt2+1:n) = Paare{z1,2}(2,schnitt2+1:n);
        end
                
    else
        Kinder{z1,1}(2,:) = Paare{z1,1}(2,:);   % Kind erhält Verteilung Chromo von Elter1, wenn PermuWhk nicht unterschritten wird
    end
    
end

