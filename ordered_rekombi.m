function [K] = ordered_crossover(K,Paare,p,n)


for z1 = 1:p
    K{z1,1} = zeros(1,n);
    
    schnitt1 = round(n*rand(1));        % zufälliger Schnitt 1
    while schnitt1 == 0
        schnitt1 = round(n*rand(1));    % Schnitt darf nicht bei Null sein
    end
    
    schnitt2 = round(n*rand(1));        % zufälliger Schnitt 2
    while schnitt1 == schnitt2 || schnitt2 == 0
        schnitt2 = round(n*rand(1));    % Schnitt darf nicht bei Schnitt1 oder Null sein
    end
    
    if schnitt2 < schnitt1
        [schnitt1,schnitt2] = deal(schnitt2,schnitt1);  % Wenn Schnitt2 kleiner als Schnitt1 ist, tauschen
    end
    
    
    K{z1,1}(1,schnitt1:schnitt2) = Paare{z1,1}(1,schnitt1:schnitt2);
    for z2 = 1:n                                % Permutationen in Elter2 hochzählen
        permu_kind = unique(K{z1,1}(1,:));      % Ermitteln welche Permutationen bereits im Kind sind
        permu_kind(permu_kind==0)=[];           % Null löschen
        
        if any(permu_kind == Paare{z1,2}(1,z2))==false         % Abfrage: Ist Permu von Elter2 schon in Kind
                            % Wenn ja (true), nächste Permu
                            % Wenn nein (false), eintragen in Kind
                                                                                
            pos_nullen = find(K{z1,1}(1,:)==0);                % Positionen der Nullen im Kinder suchen
            K{z1,1}(1,pos_nullen(1,1)) = Paare{z1,2}(1,z2);    % Wert von Elter2 in erste Null im Kind eintragen
        end
    end
end