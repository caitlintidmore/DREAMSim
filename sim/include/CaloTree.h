#ifndef CaloTree_h
#define CaloTree_h 1

#include <iostream>
#include <fstream> // for input/output files
#include <sstream> // for string stream
#include <math.h>  // for sin(x) etc.
#include <cstdlib> // for rand() on archer.
#include <iomanip> // for setw() in cout,
#include <vector>
#include <map>
#include <memory>
#include <string>

class TFile;
class TTree;
class TH1D;
class TH2D;

class CaloHit;

using namespace std;

class CaloTree
{
public:
   CaloTree(string, int argc, char **argv); // string outname
   ~CaloTree();                             // string outname
   void BeginEvent();
   void EndEvent();
   void EndJob();

   //  inout paramter handling ....
   bool setParam(string key, string val);
   float getParamF(string key);
   int getParamI(string key);
   string getParamS(string key);

   //  called fro SteppingAction...
   void accumulateHits(CaloHit aHit);
   void saveBeamXYZE(string, int, float, float, float, float);

   // for histogrming...
   std::string title;
   std::map<std::string, TH1D *> histo1D;
   std::map<std::string, TH1D *>::iterator histo1Diter;
   std::map<std::string, TH2D *> histo2D;
   std::map<std::string, TH2D *>::iterator histo2Diter;

private:
   // private functions.
   void readMacFile(string);
   vector<string> parse_line(string line);
   // std::map<std::string, std::string> mcParams;  //  MC run time parameters.
   map<string, string> mcParams; //  MC run time parameters.

   //
   void clearCaloTree();
   void analyze();

   map<int, double> make2Dhits(map<int, double> hits);
   void defineCSV(string type);
   void writeCSV(string type, map<int, double> &hits);

   string beamType;
   int beamID;
   float beamX, beamY, beamZ, beamE;
   //
   string runConfig;
   int runNumber;

   int eventCounts;
   int eventCountsALL;

   bool createNtuple;

   // hit data in csv file
   map<string, int> csvEvents; // number of events to be written to csv file.
   map<std::string, std::unique_ptr<std::ofstream>> fcsv;
   // map<string,*ofstream> csvfiles;

   // ofstream hit2DFile;  //
   // ofstream hit3DFile;  //

   // ntuple file definition...
   TFile *fout;
   TTree *tree;

   //  accumulated energyr of photons
   //  in rods
   map<int, double> rtHits; // T-slice  (nominal 50 ps/slicen), edep
   map<int, double> rzHits; // Z-slice  (nominal 2 cm/slice)  , edep
   map<int, double> rzEdep; // Z-slice  (nominal 2 cm/slice)  , edep

   // in scit fibers
   map<int, double> stHits; // T-slice  (nominal 50 ps/slicen), edep-birk
   map<int, double> szHits; // Z-slice  (nominal 2 cm/slice)  , edep-birk
   map<int, double> szEdep; // Z-slice  (nominal 2 cm/slice)  , edep

   // in sherenkov fibers
   map<int, double> ctHits; // T-slice  (nominal 50 ps/slicen), n-photons
   map<int, double> czHits; // Z-slice  (nominal 2 cm/slice)  , n-photons
   map<int, double> czEdep; // Z-slice  (nominal 2 cm/slice)  , edep

   int mRun;
   int mEvent;

   int m_run;
   int m_event;

   float m_beamMinE;  // GeV
   float m_beamMaxE;  // GeV
   float m_gridSizeX; // mm
   float m_gridSizeY; // mm
   float m_gridSizeT; // ps
   float m_caloRotationX;
   float m_caloRotationY;

   float m_calibCore;  // ph to GeV
   float m_calibHelix; // ph to GeV

   float m_calibSen;
   float m_calibSph;
   float m_calibCen;
   float m_calibCph;

   int m_beamID; //  pdg ID
   float m_beamX;
   float m_beamY;
   float m_beamZ;
   float m_beamE;
   string m_beamType;

   int m_nhits3dSS;
   vector<int> m_id3dSS;   // channel ID   xxxyyyzzz
   vector<int> m_tkey3dSS; // T-slice
   vector<int> m_zkey3dSS; // Z-slice
   vector<int> m_type3dSS;
   vector<int> m_area3dSS;
   vector<int> m_ix3dSS;
   vector<int> m_iy3dSS;
   vector<int> m_ixx3dSS;
   vector<int> m_iyy3dSS;
   vector<int> m_zslice3dSS;
   vector<int> m_tslice3dSS;
   vector<float> m_ph3dSS; // number of photons
   float m_sum3dSS;

   int m_nhits3dCC;
   vector<int> m_id3dCC; //  chnanel ID  xxxyyyttt
   vector<int> m_tkey3dCC;
   vector<int> m_zkey3dCC;
   vector<int> m_type3dCC;
   vector<int> m_area3dCC;
   vector<int> m_ix3dCC;
   vector<int> m_iy3dCC;
   vector<int> m_ixx3dCC;
   vector<int> m_iyy3dCC;
   vector<int> m_zslice3dCC;
   vector<int> m_tslice3dCC;
   vector<float> m_ph3dCC; //  number of photons
   float m_sum3dCC;
};

#endif
