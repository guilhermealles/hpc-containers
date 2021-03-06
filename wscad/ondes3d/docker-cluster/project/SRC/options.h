#ifndef OPTIONS_H
#define OPTIONS_H
#include "struct.h"

/* options generated by ./auto.sh */
static const enum typeAnelas ANLmethod=ELASTIC;
static const enum typeInterface interface=USUAL;
static const enum typeSurface surface=FREE;
static const enum typePML ABCmethod=CPML;
static const char PRMFILE[50]="./ESSAI-XML/essai.prm";
/* fixed options */
static const enum typeModel model=LAYER;
static const enum typeSource source=HISTFILE;
static const enum typeSea sea=SEA;
static const enum typeSnapshot snapType=OVELO;

static const int STATION_STEP= 10;
static const int SURFACE_STEP= 10;

#endif
