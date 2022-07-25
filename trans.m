%% Transformation einige Konfigs

%% Parameter der Fitnessfunktion
anz_fit = anz_eig + 4;  % 4 = Taktzeit, Auslastung, Pausen und Werkzeuge werden immer berechnet

% Skalieren der Gewichtung auf [0,1]
summe_gewichtung = w(1,1) + w(2,1) + w(5,1) + w(6,1) + w(7,1) + w(8,1) + w(9,1);
w(1,1) = w(1,1) / summe_gewichtung;
w(2,1) = w(2,1) / summe_gewichtung;

w(5,1) = w(5,1) / summe_gewichtung;
w(6,1) = w(6,1) / summe_gewichtung;
w(7,1) = w(7,1) / summe_gewichtung;
w(9,1) = w(9,1) / summe_gewichtung;
w(8,1) = w(8,1) / summe_gewichtung;


%% Abbruchbedingungen
time_max = abbruch_time * 60;


%% Parameter des genetischen Algorithmus
% Ersetzungsstrategie
p_el = mu_el * p;          % Anteil der Eltern für die Umweltselektion
p_ki = mu_ki * p;          % Anteil der Kinder für die Umweltselektion

% Eltern Ersetzung durch mutierte Eltern
anz_mut_eltern = round(anteil_mut_eltern*p);  






