# Masterarbeit Implementierungen

Dieses Projekt implementiert die Single Key Architektur und die MPC in TFHE innerhalb der OpenFHE für die Masterarbeit.

## Voraussetzungen

- Aktuelle jsoncpp-Bibliothek
- Aktuelle openFHE-Bibliothek

Achtung: Für das erfolgreiche Kompilieren und Erfassen der Daten für die Zeiterfassung muss die Bibliothek mit der hier beigefügten Datei `ckksrns-fhe.cpp` angepasst werden.

Diese muss in das folgende Verzeichnis eingefügt werden: 
`src/pke/lib/scheme/ckksrns/ckksrns-fhe.cpp`  

Anschließend muss die Bibliothek neu kompiliert werden.

## Nutzung

Innerhalb des Projektordners `MasterarbeitImplementierung` sind folgende Schritte notwendig:

1. Erstellen Sie einen `build`-Ordner.
2. Wechseln Sie in den `build`-Ordner.
3. Führen Sie `cmake ..` aus.
4. Führen Sie `make` aus.

Daraus erhält man die vier ausführbaren Dateien:
- BenchmarkUC1 
- BenchmarkUC1hybrid 
- BenchmarkUC2  
- BenchmarkUC2hybrid

## Source Code 

Im Verzeichnis `src` sollten `BenchmarkUC1.cpp`, `BenchmarkUC1hybrid.cpp`, `BenchmarkUC2.cpp` und `BenchmarkUC2hybrid.cpp` zur Verfügung stehen.

## Daten 

Die Daten sind synthetisch hergestellt und sind in der Datei `FleetData.json` enthalten.

# Benchmark der openFHE Bibliothek

Starten Sie das Benchmark beispielsweise mit `sudo ./BenchmarkUC1 <filename.json> <repetitions>`.

Das Benchmark variiert die Eingangsparameter:
* Rescaling Mode (FIXEDAUTO, FLEXIBLEAUTO, FLEXIBLEAUTOEXT)
* Batchsize (8, 16, 32, 64, 128 bit)
* ScaleModSize (49, 54, 59)
* FirstModSize (50, 55, 60)
* Ring Dimension (1024, 2048, 4096, 8192, 16384)

Hinweis: Die Parameter können je nach Bedarf verändert werden. In Use-case 2 wird das Security Level gesetzt und die Ringdimension automatisch ermittelt. Zeiterfassungsbefehle im Code sollten überprüft werden und können nach Belieben positioniert werden.

## Datenblatt

Die einzelnen Schritte werden gemessen und anschließend in die angegebene `.json`-Datei gespeichert.
