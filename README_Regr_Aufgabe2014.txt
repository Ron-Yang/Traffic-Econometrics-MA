Hallo Frau Solbrig,

interessanterweise haben Sie einen aehnlichen statistischen Sachverhalt wie ich bei einer Auswertung einer Umfrage der Alpenvereinssektion SBB, bei der ich eherenamtlicher "Chef-Statistiker" bin. Generell gibt es bei "hoeherwertigeren" Skalierungen auch andere Methoden, lineare Abhaengigkeiten (nicht im Sinne von Ursache und Wirkung!) zu untersuchen:

* Pearson'scher (normaler) Korrelationskoeffizient r(x,y) bei zwei kardinalskalierten Groessen x und y

* Spearmann'sche Rangkorrelation rs(x,y) bei zwei ordinalskalierten oder einer ordinal- und einer kardinalskalierten Groesse. Es gilt

rs(x,y)=r(R(x), R(y))

In Worten: Der Rangkorrelationskoeffizient von x und y ist gleich dem Pearsonschen Korrelationskoeffizient der Rangziffern R(x) und R(y) von x bzw. y (der niedrigste Wert hat Rangziffer 1, der hoechste die Ziffer n oder umgekehrt).


Beide Korrelationen kann man sowohl direkt (deskriptiv) als Mass der linearen Abhaengigkeit anwenden, mit dem bekannten Wertebereich r(x,y) in [-1,1]. Man kann aber auch  induktiv einen Korrelationstest darauf aufbauen, also den Test  der Nullhypothese "es gibt keine lineare Abhaengigkeit, r=0" ( siehe mein Oekonometrie- oder Statistik-Skript): Die Testvariable

T=r(x,y)*sqrt(n-2)/sqrt(1-r(x,y)) sim T(n-2)

ist (bei Zutreffen der Nullhypothese und gaussverteilter Fehlerterme) student-t verteilt mit n-2 Freiheitsgraden

Wichtig:

(1) Diese Tests sind nur bei vermuteter *linearer* Abhaengigkeit sinnvoll. Hat man bspw eine parabelfoermige Abhaengigkeit (grosse Werte von y fur kleine und grosse x, kleine fuer mittlere x), werden die Korrelationstests insignifikant, auch wenn es eine hochsignifikante, aber eben *nichtlineare* Abhaengigkeit gibt. Der Unabhaengigkeitstest findet diese aber immer, egal wie die Skalierung ist!

(2) Im Vierfelderfall (sowohl x als auch y haben nur je zwei Auspraegungen) sind alle Tests aequivalent

(3) Auch fuer den Unabhaengigkeitstest gibt es eine deskriptiv anwendbare Zahl: Der

Determinationskoeffizient r_det^2=Q/(n*min(Kx-1, Ky-1)) in [0,1]

mit Q der Testvariable des Chi^2 Tests, n=Datenumfang und Kx und Ky die Klassenzahlen von x bzw. y

Im Vierfelderfall gilt r_det^2=r^2

Also: Wenn Sie vermuten, dass eine lineare Abhaengigkeit dominant ist, z.B. dass Leute umso positiver bewerten, je laenger  die Reiseweite ist, nehmen Sie den Korrelationstest (bzgl. der Rangziffern), ansonsten (z.B. Leute mit kurzer und langer Reiseweite bewerten positiv, welche mit mittlerer negativ) den Unabhaengigkeitstest. Fuer diesen gibt es natuergemaess (ausser r_det^2) keine direkte Analysezahl, sondern Sie muessen bei Signifikanz direkt die Kreuztabellen diskutieren

Wie Sie bei so wenigen Auspraegungen allerdings bei n=609 auf Kreuztabellenelemente mit h(x_i,y_j)<5 (Unabhaengigkeitstest nicht mehr anwendbar) kommen, ist mir raetselhaft. Fassen Sie im Zweifelsfall Auspraegungen zusammen (z.B. auf drei Entfernungsklassen und vier Bewertungsklassen; falls dann immer noch zu wenige "ganz schlecht" bewerten, fassen Sie ausserdem die beiden schlechtesten Bewertungsklassen zusammen): Weniger ist manchmal mehr!

Schliesslich koennen Sie einen Zusammenhang zwischen einer kardinalskalierten und einer ordinal- oder nominalskalierten Groesse auch mit der Regressionsrechnung untersuchen. Vorteil: Sie verschenken nicht unnoetig Genauigkeit durch das sonst notwendige Kategorisieren der kardinalskalierten Variablen in Klassen. Wichtig dabei ist, dass die endogene Variable immer die kardinalskalierte sein muss, auch wenn die vermutete Ursache-Wirkungs-Richtung andersherum ist, also z.B hier:

Endogene Variable y=Reiseweite in R+,
Exogene Variable x=Bewertung in {1,2,3,4}
Regressionsfunktion haty(x)=beta0+beta1*delta(x,1)+beta2*delta(x,2)+beta3*delta(x,3)

mit delta(x,j)=1 falls x=j und 0 sonst (x=4 ist Referenzauspraegung, siehe Skript).

Dies funzt auch bei beliebiger Nichtlinearitaet. Ist bspw beta1>0 sowie >beta2 und >beta3, so reisen Leute, die mit '1' bewerten, bevorzugt lang. Sind hingegen beta1, beta2, beta3 alle <0, so sind die Langreisenden bevorzugt die, welche mit '4' bewerten.

Hingegen ist der Test auf Null des Parameters beta1 der *linearen* Funktion haty(x)=beta0+beta1*x identisch zum Test auf eine verschwindende Korrelation, r(x,y)=0

Viele Gruesse,

Martin Treiber


On 01/07/2014 03:16 PM, Anka Solbrig wrote:
> Sehr geehrter Dr. Treiber,
>
> vielen Dank für Ihre hilfreichen Ausführungen.
>
> Noch eine Frage, ob mein Vorgehen richtig ist.
> Die Chemnitzer Verkehrs-AG hat 15.000 Abo-Kunden, davon wurden 609 telefonisch befragt und es liegen Daten vor.
> Unter anderem wurden Bewertungen nach Noten (1-4, also ordinal) z. B. für Takt, Linien- und Streckennetz usw. abgefragt. Aus weiteren Antworten wurden die Reiseweiten und Reisezeiten (metrisch) usw. berechnet.
> Für jeden der 609 befragten Abo-Kunden gibt es also die Bewertungen und mittleren Reiseweiten.
> Ich möchte prüfen, ob es einen Zusammenhang zwischen der Bewertung und z. B. der Reiseweite gibt.
> Also ob ein Zusammenhang zwischen Bewertung des Taktes und der mittleren Reiseweite gibt,
>
> Ist der Chi-Quadrat-Test auch hier der richtige? Oder gibt es einen passenderen für das Niveau von ordinal/metrischen Daten?
> Den Chi-Quadrat-Test habe ich bereits bei einigen Versuchen angewendet. Häufigstes Problem ist wieder, dass die erwarteten Häufigkeiten zu oft zu klein sind.
> Kann ein Ergebnis auch sein, dass der Zusammenhang nur schwer nachgewiesen werden kann sich aber vermuten lässt?
> Ich weiß nicht, welche Ausmaße dasTesten annehmen soll.
>
> Besten Dank und freundliche Grüße
>
> Anka Solbrig
>
>
>
>
> Am 19. Dezember 2013 13:19 schrieb Martin Treiber <treiber@vwi.tu-dresden.de>:
>
>     Hallo Frau Solbrig,
>
>     es ist richtig, dass bei erwarteten Haeufigkeiten unter 5 der chi^2-Test nicht genau genug ist, da dann der Unterschied seiner Naeherungen zum exakten Abzaehlen mittels hypergeometrischer Verteilung (=exakter Fisher) zu gross ist. Meine eigenen Tests in Faellen, wo man beide berechnen kann, besagen, dass der chi^2 Test sich dann zu weit aus dem Fenster lehnt, d.h. kleinere p-Werte rausgibt. Dies gilt fuer den Vierfelderfall sogar, wenn alle erw Haeufigkeiten > 5 sind. Andererseits wird aber auch gesagt, dass der exakte Test wegen der Granularitaet bei kleinen absoluten Haeufigkeiten zu konservativ ist (zu grosse p-Werte ergibt). Dennoch sollte man wohl immer den exakten Test nehmen, wenn die erwarteten Haeufigkeiten kleiner 5 sind oder es nur vier Felder gibt.
>
>     Konkret zu Ihren Fragen:
>
>     * der exakte Test ist fuer beliebige n X m Kreuztabellen durch elementarer Kombinatorik definiert. Allerdings gibt es bei Kreuztabellen mit vielen Zeilen und/oder Spalten so viele Faelle zu betrachten, dass man die Kombinatorik nicht mehr praktikabel ausrechnen kann, deshalb die langen Berechnungszeiten in Ihrem Programm.
>
>     * Falls die Berechnung zu lange dauert, sollte man Zellen zusammenfassen (z.B. mehrere Stadteile zu einem Gebiet). Dann ist einerseits die exakte Rechnung schneller und man kann anderrseits wieder den Chi^2-Test anwenden. Dessen ungeachtet sehe ich keinen Sinn, sehr kleine Teil-Stichproben (z.B. 10 Leute aus Stadteil x) herauszufiltern. Selbst wenn man dann nur 2 Frauen und 8 Maenner hat, ist aufgrund der kleinen Stichprobe nicht zu folgern, dass Maenner systematisch bevorzugt ausgewaehlt wurden.
>
>     Viele Grüße
>
>      Martin Treiber
>
>     ----------------------------------------------------------
>     Dr. Martin Treiber
>              Institute for Transport & Economics, TU Dresden
>     	 Chair of Traffic Modelling, Econometrics, and Statistics
>              Falkenbrunnen, Zimmer 104 (two floors up from the entrance!)
>              Würzburger Str. 35, D-01062 Dresden
>              treiber@vwi.tu-dresden.de,
>              www.mtreiber.de
>              phone/fax:     +49 (351) 463 36794 / 36809
>      ---------------------------------------------------------
>
>     On 12/18/2013 02:19 PM, Anka Solbrig wrote:
>>     Sehr geehrter Dr. Treiber,
>>     vielen Dank für Ihre schnelle und vor allem auch hilfreiche Antwort.
>>     Sie hat mir im ersten Moment sehr geholfen.
>>     Jedoch bin ich heute bei der Anwendung des Chi-Quadrat-Tests auf das Problem gestoßen, dass die erwarteten Häufigkeiten bei fast allen Untersuchungen in mehr als 20 % der Fälle kleiner als 5 sind. Ala Alternative habe ich in der Literatur den exakter Fisher-Test gefunden. Hier sehe ich allerdings das Problem, dass meiner Meinung nach keine 2x2 Kreuztabelle vorliegt.
>>     Zur Verdeutlichung. Ich möchte den Zusammenhang zwischen soziodemografischen Merkmalen untersuchen. Gibt es beispielsweise eine Abhängigkeit zwischen den Stadtteilen von Chemnitz (39 Möglichkeiten) und den Geschlechtern der Einwohner (2 Möglichkeiten). Aufgrund der Vielzahl der Felder liegen die erwarteten Häufigkeiten in knapp 50 % der Fälle unter 5. Ich bin der Meinung den Fisher-Test kann man nicht anwenden, bzw. ist dies in SPSS mit einer sehr langen Berechnungszeit verbunden, wo der Arbeitsspeicher des Rechners nicht ausreicht. Welcher Test ist der richtige?
>>     Nochmals besten Dank.
>>
>>     Viele Grüße
>>
>>     Anka Solbrig
>>
>>     Am 13. Dezember 2013 15:18 schrieb Martin Treiber <treiber@vwi.tu-dresden.de>:
>>
>>         Sehr geehrte Frau Solbrig,
>>
>>         der chi^2-Test ist genau der Richtige! Wichtig dabei ist, dass Sie ggf. auch nicht-nominalskalierte Merkmale (Bspw. das Alter) in Kategorien bzw. Klassen einteilen.
>>
>>         Zur Frage Stichprobe oder Grundgesamtheit: Zunaechst koennen Sie allein aus der Stichprobe den chi^2-Test genau so anwenden, wie er im Statistik-Skript steht. Haben Sie M Merkmale und wollen Sie die Abhaengigkeit jedes Merkmals von jedem anderen untersuchen, ergibt das M*(M-1)/2 jeweils getrennte chi^2-Tests. Dabei sollten Sie durchaus hohe Ansprueche an die Signifikanz stellen (z.B. alpha bzw. p<0.1 %), um statistische Selektionsfehler zu vermeiden: Bei 20 chi^2-Tests wuerde  selbst bei idealer Unabhaengigkeit im Mittel einer der Tests auf dem 5%-Niveau "anschlagen".
>>
>>         Sie koennen aber auch, fuer die Merkmale die in der Grundgesamtheit vorliegen,  die Stichprobe mit der Grundgesamtheit vergleichen und sie so (bzgl. des jeweiligen Merkmals) auf Repraesentativitaet testen. Dies ist der chi^2-Homogenitaetstest, der nichts anderes ist als ein chi^2-Unabhaengigkeitstest mit der zweitem Variablen gleich dem Index der Datenquelle (z.B. 1 fuer die Grundgesamtheit, 2 fuer die Stichprobe). Um valide Aussagen aus der Stichprobe zu erhalten, sollte der Homogenitaetstest moeglichst *nicht* anschlagen, also p-Werte oberhalb von 5% liefern (oder genauer: Der Anteil der Homogenitaetstests mit p-Werten unterhalb 5% sollte 5% nicht deutlich ueberschreiten).
>>
>>         Wichtig: Der Unabhaengigkeitstest sagt nichts ueber die Richtung der Abhaengigkeit aus. Dazu ist ggf (also wenn ein Test "anschlaegt") die direkte Analyse der jeweiligen Kreuztabelle noetig. Ausserdem wichtig: Der Unabhaengigkeitstest sagt nichts ueber Ursache-Wirkung aus. Diese kann von X nach Y, von Y nach X gehen, oder die Ursache ist eine dritte Variable.
>>
>>         Ich hoffe, das hilft Ihnen weiter.
>>
>>         Gruss,
>>
>>         Martin Treiber
>>
>>         ----------------------------------------------------------
>>         Dr. Martin Treiber
>>                  Institute for Transport & Economics, TU Dresden
>>                  Chair of Traffic Modelling, Econometrics, and Statistics
>>                  Falkenbrunnen, Zimmer 104 (two floors up from the entrance!)
>>                  Würzburger Str. 35, D-01062 Dresden
>>                  treiber@vwi.tu-dresden.de,
>>                  www.mtreiber.de
>>                  phone/fax:     +49 (351) 463 36794 / 36809
>>          ---------------------------------------------------------
>>
>>         On 12/13/2013 02:40 PM, Anka Solbrig wrote:
>>
>>             Sehr geehrte Dr. Treiber,
>>
>>             ich schreibe derzeit meine Masterarbeit bei der Chemnitzer Verkehrs-AG. Im Bachelorstudium besuchte ich bei Ihnen die Statistik-Vorlesungen. Hauptbestandteil ist die Auswertung von Daten zu den Stammkunden bei der CVAG.
>>             Zum einen gibt es Informationen über alle 15.000 Abokunden. Hauptsächlich arbeite ich aber mit Daten zu 600 Personen, die im Rahmen einer telefonischen Befragung interviewt worden. Zu diesen 600 Personen stehen soziodemografische Angaben zur Verfügung, z. B. der Stadtteil, aus dem sie stammen oder, ob ein PKW zur Verfügung steht. Die meisten Merkmale sind nominal skaliert. Jetzt möchte ich prüfen, ob Abhängigkeiten zwischen den Merkmalen bestehen. Z. B. wäre es interessant zu wissen, ob der Stadtteil mit der PKW-Verfügbarkeit zusammenhängt.
>>
>>             Meine Frage ist nun, ob für diese Problemstellung der Chi Quadrat Test anwendbar ist.Das Datenniveau ist meiner Recherche zufolge ausreichend. Jedoch bin ich nicht sicher, ob mehr Informationen über zu erwartende Häufigkeiten in der Grundgesamtheit bekannt sein müssten. In diesem Zusammenhang stellt sich die Frage, ob die Grundgesamtheit in dem Falle die 600 befragten Personen oder alle 15.000 Stammkunden sind.
>>
>>             Kann ich den Chi Quadrat Test bei der Problematik anwenden? Wenn nicht, können Sie einen anderen Test empfehlen?
>>
>>             Für Ihre Antwort besten Dank und eine schöne Weihnachtszeit.
>>
>>             Mit freundlichen Grüßen
>>
>>             Anka Solbrig 
>>

