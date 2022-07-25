function [Eltern,Population,Pointer,SUS_Fitness] = SUS_elternselektion (Population,p,fitness_typ)

%%
% for z1 = 1:p
%     Population{z1,9} = z1;
% end

P_sort = sortrows(Population,-fitness_typ);    % P nach Fitness sortieren

GesamtFitness = 0;          % Gesamte Fitness bestimmen
for z1 = 1:p
    GesamtFitness = GesamtFitness + P_sort{z1,fitness_typ}(1,1);
end


Pointer_Range = GesamtFitness / (2*p);  % Abstand der Pointer zueinander
start_sus = Pointer_Range * rand(1);    % Start der Pointer


Pointer = zeros(2*p,1);     % Position der Pointer berechnen
for elt=1:2*p
    Pointer(elt,1) = start_sus + (elt-1)*Pointer_Range; % -1, da start_sus schon der erste pointer ist
end


SUS_Fitness = zeros(p,1);

SUS_Fitness (1,1) = P_sort{1,fitness_typ}(1,1);   % aufsummieren der Fitness für Fitness Skala
for z1 = 2:p
    SUS_Fitness(z1,1) = P_sort{z1,fitness_typ}(1,1) + SUS_Fitness(z1-1,1);
end

Eltern = cell(p,2);

for elt = 1:2*p     % Pointer hochzählen
    for z1 = 1:p    
        
        if z1 == 1
            Eltern{elt,1} = P_sort{z1,1};   % Eigentlich Bedingung notw.
            Eltern{elt,2} = z1;             % z1 wird immer hinzugefügt
        end
        
        if z1 > 1
            if Pointer(elt,1) <= SUS_Fitness(z1,1) && Pointer(elt,1) > SUS_Fitness(z1-1,1)
                Eltern{elt,1} = P_sort{z1,1};   % wenn hier allerdings zutreffend, wird z1 überschrieben
                Eltern{elt,2} = z1;
            end
        end 
    end
end
