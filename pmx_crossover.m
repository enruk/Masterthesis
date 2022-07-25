function [Kinder] = pmx_crossover(Kinder,Paare,p,n,PermuCrossWhk)

%%
for z1=1:p
    
    random = rand(1);
    if random < PermuCrossWhk
        Kinder{z1,1} = zeros(1,n);

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

        Kinder{z1,1}(1,schnitt1:schnitt2) = Paare{z1,1}(1,schnitt1:schnitt2);   % Abschnitt eintragen in das Kind


        for position = schnitt1:1:schnitt2              % Abschnitt hochzählen
            place = position;                           % Place zeigt die Spalte an

            if any(Paare{z1,1}(1,schnitt1:schnitt2) == Paare{z1,2}(1,position)) == 0
                while place >= schnitt1 && place <= schnitt2
                    platzhalter1 = Paare{z1,1}(1,place);
                    place = find(Paare{z1,2}(1,:) == platzhalter1);
                end
                Kinder{z1,1}(1,place)=Paare{z1,2}(1,position);
            end

        end


        % Letzte Nullen nach Reihenfolge der Permu in Elter2 ersetzen
        for z2 = 1:n                                % Permutationen in Elter2 hochzählen
        permu_kind = unique(Kinder{z1,1}(1,:));         % Ermitteln welche Permutationen bereits im Kind sind
        permu_kind(permu_kind==0)=[];                   % Null löschen

            if any(permu_kind == Paare{z1,2}(1,z2))==false         % Abfrage: Ist Permu von Elter2 schon in Kind
                                % Wenn ja (true), nächste Permu
                                % Wenn nein (false), eintragen in Kind

                pos_nullen = find(Kinder{z1,1}(1,:)==0);                % Positionen der Nullen im Kinder suchen
                Kinder{z1,1}(1,pos_nullen(1,1)) = Paare{z1,2}(1,z2);    % Wert von Elter2 in erste Null im Kind eintragen
            end
        end
        
        
    else
        Kinder{z1,1}(1,:) = Paare{z1,1}(1,:);   % Kind erhält Permu Chromo von Elter1, wenn PermuWhk nicht unterschritten wird
    end
    
end



