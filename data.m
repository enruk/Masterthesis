[xlsfile,path2xls] = uigetfile('*.xls','Bitte Datei auswählen');
fid = fopen (fullfile(path2xls,xlsfile),'r'); 
[num,txt,raw] = xlsread(xlsfile);

n = max(num(:,1));          % n: Anzahl der Operationen
T = zeros(n,2);             % T Matrix der Prozesszeiten
WZ = zeros(n,1);             % T Matrix der Prozesszeiten
E = zeros(n,2);             % E Matrix des 1. Eignungsgrad
E2 = zeros(n,2);            % E Matrix des 2. Eignungsgrad
E3 = zeros(n,2);            % E Matrix des 3. Eignungsgrad
E4 = zeros(n,2);            % E Matrix des 4. Eignungsgrad
E5 = zeros(n,2);            % E Matrix des 5. Eignungsgrad

logicalpos = strcmp(raw, 'Werkzeug');
[~,sNr] = find(logicalpos);
posWZ = sNr;        % sNr: SpaltenNummer

logicalpos = strcmp(raw, 'Dauer M');
[~,sNr] = find(logicalpos);
posTM = sNr;

logicalpos = strcmp(raw, 'Dauer R');
[~,sNr] = find(logicalpos);
posTR = sNr;

logicalpos = strcmp(raw, 'KM');
[~,sNr] = find(logicalpos);
posKM = sNr;

logicalpos = strcmp(raw, 'KR');
[~,sNr] = find(logicalpos);
posKR = sNr;

logicalpos = strcmp(raw, 'KM2');
[~,sNr] = find(logicalpos);
posKM2 = sNr;

logicalpos = strcmp(raw, 'KR2');
[~,sNr] = find(logicalpos);
posKR2 = sNr;

logicalpos = strcmp(raw, 'KM3');
[~,sNr] = find(logicalpos);
posKM3 = sNr;

logicalpos = strcmp(raw, 'KR3');
[~,sNr] = find(logicalpos);
posKR3 = sNr;

logicalpos = strcmp(raw, 'KM4');
[~,sNr] = find(logicalpos);
posKM4 = sNr;

logicalpos = strcmp(raw, 'KR4');
[~,sNr] = find(logicalpos);
posKR4 = sNr;

logicalpos = strcmp(raw, 'KM5');
[~,sNr] = find(logicalpos);
posKM5 = sNr;

logicalpos = strcmp(raw, 'KR5');
[~,sNr] = find(logicalpos);
posKR5 = sNr;


% Matrix der Werkzeuge
WZ(:,1) = num(:,posWZ);

% befüllen der Matrizen für Zeit und Eignungsgrad
T(:,1) = num(:,posTM);
T(:,2) = num(:,posTR);

if anz_eig >= 1
    E(:,1) = num(:,posKM);
    E(:,2) = num(:,posKR);
end

if anz_eig >= 2
    E2(:,1) = num(:,posKM2);
    E2(:,2) = num(:,posKR2);
end

if anz_eig >= 3
    E3(:,1) = num(:,posKM3);
    E3(:,2) = num(:,posKR3);
end

if anz_eig >= 4
    E4(:,1) = num(:,posKM4);
    E4(:,2) = num(:,posKR4);
end

if anz_eig >= 5
    E5(:,1) = num(:,posKM5);
    E5(:,2) = num(:,posKR5);
end

% Generieren der Präzedenzmatrix

Prae = zeros(n,n);          % Präzedenzmatrix (oder auch Vorrangmatrix genannt)
A_ges = zeros(n,n);         % Matrix der Vorgänger jeder Operation       

        % wichtige Spalten bestimmen und raw Matrix anpassen
logicalpos = strcmp(raw, 'Nr');     % suche Spalte Nr
[~,sNr] = find(logicalpos);           
k = 1;
while raw{k,sNr}(1,1)~= 1           % suche '1' in Spalte Nr
    k = k + 1;
end

logicalpos = strcmp(raw, 'Vor');    % suche Spalte Vor
[~,sVor] = find(logicalpos); 

raw2 = raw(k:end,1:end);            % kürzen der raw Matrix sodass Spalte 'Nr' mit 1 beginnt
sz = size(raw2);
zmax = sz(1,1);

        % Füllen der Präzidenzmatrix Prae
for z = 1:1:zmax
    str=raw2{z,sVor};   
    if isnumeric(str) == 0          % aussortieren der NaN Zeilen (keine Präzdenz beim betreffenden Prozess)
        prozessnr = z;
        index = str2num(regexprep(raw2{z,sVor}, '\D+', ' '));   % Zahlen aus str bestimmen
        sz = size(index);
        N = sz(1,2);    
                                    % Zahlen aus dem Index Vektor in
                                    % Schleife abgehen
        for k = 1:N
            vorprozessnr = index(1,k);
            Prae(prozessnr,vorprozessnr) = 1;                   % Prae beschreiben
            A_ges(z,k) = index(1,k);
        end
    end
end


%%

clearvars raw raw2 sNr str txt xlsfile num vorprozessnr zmax sz sVor prozessnr path2xls N logicalpos
clearvars posKM posKM2 posKM3 posKM4 posKM5 posKR posKR2 posKR3 posKR4 posKR5 posTM posTR

fclose('all');

