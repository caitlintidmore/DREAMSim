#ifndef CaloHit_h
#define CaloHit_h 1

#include <vector>

#include "CaloID.h"

class CaloHit
{
public:
   CaloHit();
   ~CaloHit();

   void print();

   CaloID caloid;
   int pid;
   int trackid;
   double x;
   double y;
   double z;
   double steplength;
   double globaltime;
   double edep;
   double edepbirk;
   double ncer;    // number of cerenkov photons
   double ncercap; // number of cerenkov photons
};

#endif
