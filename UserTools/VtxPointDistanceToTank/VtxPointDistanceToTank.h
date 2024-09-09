#ifndef VTXPOINTDISTANCETOTANK_H
#define VTXPOINTDISTANCETOTANK_H



#include <string>
#include <iostream>
#include "Tool.h"


#include "TApplication.h"
#include <Math/PxPyPzE4D.h>
#include <Math/LorentzVector.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "ADCPulse.h"
#include "Waveform.h"
#include "CalibratedADCWaveform.h"
#include "Hit.h"
#include "RecoDigit.h"
#include "ANNIEalgorithms.h"
#include "TimeClass.h"
#include "BeamStatus.h"
#include "ANNIEGeometry.h"	
#include "Position.h"

/**
 * \class VtxPointDistanceToTank
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: Christian nguyen s $
* $Date: 2024 $
* Contact: cnguyen@fnal.gov
*/
class VtxPointDistanceToTank: public Tool {


 public:



  VtxPointDistanceToTank(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  bool GetSimpleRECO_Vertex();
  bool GetTRUE_Vertex();

 private:

void ResetVariables();
  //General variables
  bool isData;
  bool hasGenie;
  bool makeSimple;
//Geometry *geom = nullptr;

  Position SimpleRecoVtx;

  double fTrueVtxX;
  double fTrueVtxY;
  double fTrueVtxZ;
  double fTrue_DistanceToEdge;

  double fRecoVtxX;
  double fRecoVtxY;
  double fRecoVtxZ;
  double fReco_DistanceToEdge;

  double fSimpleVtxX;
  double fSimpleVtxY;
  double fSimpleVtxZ;
  double fSimpleVtx_DistanceToEdge;

  // Trigger-level information
  std::map<std::string,bool> fDataStreams;
  int fTriggerword;
  int fTankMRDCoinc;
  int fNoVeto;
  int fHasTank;
  int fHasMRD;

  /// \brief Integer that determines the level of logging to perform
  int verbosity = 0;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;		
  int get_ok;	

};


#endif
