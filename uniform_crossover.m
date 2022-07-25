function [Kinder]=uniform_crossover(Kinder,Paare,p,n,co_unif_whk,AlloCrossWhk)

% Unterscheide zur Uniform Mutation bei bin�ren Repr�sentationen
    % Treffen zwei Indi mit �hnlichen Genen aufeinander, tritt keine
    % gro�e Ver�nderung ein: Gerade zum Ende sinkt so die Diversit�t

for z1=1:p
    
    random = rand(1);
    if random < AlloCrossWhk
        
        Kinder{z1,1}(2,:) = zeros(1,n);     % Zeile 2 beim Kind erstmal mit nullen f�llen

        for z2=1:n
            random = rand(1);               % Zufallszahl bestimmen
            if random > co_unif_whk
                Kinder{z1,1}(2,z2) = Paare{z1,1}(2,z2); % wenn random gr��er ist, dann Gen von Elter1
            elseif random <= co_unif_whk
                Kinder{z1,1}(2,z2) = Paare{z1,2}(2,z2); % wenn random kleiner gleich ist, dann Gen von Elter2
            end
        end

    else
        Kinder{z1,1}(2,:) = Paare{z1,1}(2,:);   % Kind erh�lt Verteilung Chromo von Elter1, wenn PermuWhk nicht unterschritten wird
    end
    
end



























