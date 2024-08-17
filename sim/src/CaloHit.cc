#include "CaloHit.h"

#include <iostream> // for cout

CaloHit::CaloHit()
{
   x = 0.0;
}

CaloHit::~CaloHit()
{
}

// ------------------------------------------------------------------
void CaloHit::print()
{
   std::cout << "CaloHit pid " << pid << " xyz(" << x << "," << y << "," << z << ") L " << steplength;
   std::cout << " T " << globaltime << "  edep " << edep << "  ncer " << ncer;
   std::cout << " " << std::endl;
   caloid.print();
}
