# medicalRadar
Pythonkode for signalbehandling av medinsk radardata. Mye av dette ble skrevet rundt 2011-2012 for daværende pythonversjon, og har ikke blitt oppdatert siden. Koden er designet for å jobbe med noen spesifikke formater, og for spesifikke formål.

## Measurement.py
Målingsobjekt, inneholder funksjoner for å lese inn data, noe prossessering, og plotting av data.

## accelerometer_display.py
Kort script for å lese inn akselerometerdata, konvertere til avstander, og plotte disse.

## chirpz.py
Chirp-Z-transformkode hentet fra nett. Modifisert for å passe formålet, og for å returnere kun delene av det komplekse plan som er interessante for målingen.

## contact_CW.py
Script for å regne på kombinerte ECG, stetoskop og radardata og plotting av resultatene fra dette.

## plot_dielectrics.py
Script for å regne på dielektrisitetsmålinger og fremvisning av dette.

## processing.py
En samling digital signalbehandlingsmetoder til prosessering av data.
