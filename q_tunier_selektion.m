function PK = q_tunier_selektion (PK,p_sel,q,champ,umwelt_fitness,champion_confi)


for z1=1:p_sel
    PK{z1,9}(1,1) = 0;      % Anzahl der Siege in Spalte 9 zurücksetzen für alle
end

if champion_confi == 1      % Champion-Funktion ein
    siege_champ = 10;       % Siege des Champion hochsetzen
elseif champion_confi == 0  % Champion-Funktion aus
    siege_champ = 0;        % Siege bleiben bei 0
end

PK{champ,9}(1,1) = siege_champ;      % sicherstellen, dass Champion Tunier übersteht    


for index = 1:q         % Tunier 1-q hochzählen
    
    tuniergegner(:,1) = randperm(p_sel);
    
    for z1=1:(p_sel/2)            % Matches 1 bis p hochzählen
        gegner1 = find(tuniergegner(:,1)==2*z1-1);
        gegner2 = find(tuniergegner(:,1)==2*z1);
        
        if umwelt_fitness == 1
            % Fitnessvergleich mit normaler Fitness
            if PK{gegner1,2}(1,1) > PK{gegner2,2}(1,1)
                PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
            elseif PK{gegner1,2}(1,1) < PK{gegner2,2}(1,1)
                PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
            elseif PK{gegner1,2}(1,1) == PK{gegner2,2}(1,1)
                PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 0;
                PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 0;
            end
        end
        
        if umwelt_fitness == 2
            % Fitnessvergleich mit shared Fitness
            if PK{gegner1,6}(1,1) > PK{gegner2,6}(1,1)
                PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
            elseif PK{gegner1,6}(1,1) < PK{gegner2,6}(1,1)
                PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
            elseif PK{gegner1,6}(1,1) == PK{gegner2,6}(1,1)
                PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 0;
                PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 0;
            end
        end
        
        
        if umwelt_fitness == 3
            % Rangvergleich mit Crowding Distance
            if PK{gegner1,5}(1,1) > PK{gegner2,5}(1,1)
                PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
            elseif PK{gegner1,5}(1,1) < PK{gegner2,5}(1,1)
                PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
            elseif PK{gegner1,5}(1,1) == PK{gegner2,5}(1,1)

                hamming1 = mean(PK{gegner1,7}(:,1));
                delta_zeit1 = mean(PK{gegner1,7}(:,2));
                %crowding1 = (hamming1+delta_zeit1)/2;
                
                
                hamming2 = mean(PK{gegner2,7}(:,1));
                delta_zeit2 = mean(PK{gegner2,7}(:,2));
                %crowding2 = (hamming2+delta_zeit2)/2;

%                 if crowding1 < crowding2
%                     PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
%                 elseif crowding1 > crowding2 
%                     PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
%                 elseif crowding1 == crowding2
%                     PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
%                     PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
%                 end
%                 

                if (hamming1 <= hamming2 && delta_zeit1 < delta_zeit2) || (hamming1 < hamming2 && delta_zeit1 <= delta_zeit2)
                    PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
                    
                elseif (hamming1 >= hamming2 && delta_zeit1 > delta_zeit2) || (hamming1 > hamming2 && delta_zeit1 >= delta_zeit2) 
                    PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
                    
                elseif hamming1 == hamming2 && delta_zeit1 == delta_zeit2
                    PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
                    PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
                    
                elseif hamming1 > hamming2 && delta_zeit1 < delta_zeit2 || hamming1 < hamming2 && delta_zeit1 > delta_zeit2
                    PK{gegner1,9}(1,1) = PK{gegner1,9}(1,1) + 1;
                    PK{gegner2,9}(1,1) = PK{gegner2,9}(1,1) + 1;
                    
                end
            end
        end
        
        
    end
    % alle Matches fertig
    % dann nächstes Tunier mit neuen Gegnern
    
end