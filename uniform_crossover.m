function [Kinder]=uniform_crossover(Kinder,Paare,p,n,co_unif_whk,AlloCrossWhk)

% Unterscheide zur Uniform Mutation bei binären Repräsentationen
    % Treffen zwei Indi mit ähnlichen Genen aufeinander, tritt keine
    % große Veränderung ein: Gerade zum Ende sinkt so die Diversität

for z1=1:p
    
    random = rand(1);
    if random < AlloCrossWhk
        
        Kinder{z1,1}(2,:) = zeros(1,n);     % Zeile 2 beim Kind erstmal mit nullen füllen

        for z2=1:n
            random = rand(1);               % Zufallszahl bestimmen
            if random > co_unif_whk
                Kinder{z1,1}(2,z2) = Paare{z1,1}(2,z2); % wenn random größer ist, dann Gen von Elter1
            elseif random <= co_unif_whk
                Kinder{z1,1}(2,z2) = Paare{z1,2}(2,z2); % wenn random kleiner gleich ist, dann Gen von Elter2
            end
        end

    else
        Kinder{z1,1}(2,:) = Paare{z1,1}(2,:);   % Kind erhält Verteilung Chromo von Elter1, wenn PermuWhk nicht unterschritten wird
    end
    
end



























