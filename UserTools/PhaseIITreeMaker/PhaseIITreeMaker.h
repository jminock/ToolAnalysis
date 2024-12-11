#ifndef PhaseIITreeMaker_H
#define PhaseIITreeMaker_H

#include <string>
#include <iostream>

#include "Tool.h"
// ROOT includes
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

#include "GenieInfo.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/JamesRandom.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/Kinematics.h"
//#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/AppInit.h"
#include <Tools/Flux/GSimpleNtpFlux.h>
#include <Tools/Flux/GNuMIFlux.h>
#include <Framework/GHEP/GHepUtils.h>               // neut reaction codes
#include <Framework/ParticleData/PDGLibrary.h>
#include <Framework/ParticleData/PDGCodes.h>
#include <Framework/Ntuple/NtpMCEventRecord.h>
#include <Framework/Ntuple/NtpMCTreeHeader.h>
#include <Framework/Conventions/Constants.h>
#include <Framework/GHEP/GHepParticle.h>
#include <Framework/GHEP/GHepStatus.h>
#include <TParticlePDG.h>
#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class PhaseIITreeMaker: public Tool {


 public:

  PhaseIITreeMaker();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  /// \brief Load MCTruth information into ROOT file variables
  //
  /// Function loads the Muon MC Truth variables with the information
  /// Saved in the RecoEvent store.
  bool FillMCTruthInfo();
  bool FillTankRecoInfo();
  int LoadMRDTrackReco(int SubEventNumber);
  void LoadAllMRDHits(bool IsData);
  void FillRecoDebugInfo();
  void FillSimpleRecoInfo();
  void FillMuonFitterInfo();
  void FillRingCountingInfo();
  void FillWeightInfo();
  void FillTruthRecoDiffInfo(bool got_mc, bool got_reco);
  void LoadDigitHits();

  /// \brief Summary of Reconstructed vertex
  void RecoSummary();
  void LoadTankClusterHits(std::vector<Hit> cluster_hits);
  void LoadTankClusterHitsMC(std::vector<MCHit> cluster_hits,std::vector<unsigned long> cluster_detkeys);
  bool LoadTankClusterClassifiers(double cluster_time);
  void LoadAllTankHits(bool IsData);
  void LoadSiPMHits();

  bool is_meson_or_antimeson( int pdg_code );
  void FillCC0pi_categorize_event(bool IsMC, bool insideFV, int isCC, int TrueNuPDG, std::vector<int> Event_PDGs, std::vector<double> Event_P, int TrueFSLPdg, double inputMC_Pmu, double mcmrdtracklength, double mcangle);
  bool FidVol(double x, double y, double z);  
  
 private:

  //General variables
  bool isData;
  bool hasGenie;

  std::map<int,std::string>* AuxChannelNumToTypeMap;
  std::map<int,double> ChannelKeyToSPEMap;

  std::map<int,unsigned long> pmtid_to_channelkey;
  std::map<unsigned long, int> channelkey_to_pmtid;
  std::map<unsigned long, int> channelkey_to_mrdpmtid;
  std::map<int, unsigned long> mrdpmtid_to_channelkey_data;
  std::map<unsigned long, int> channelkey_to_faccpmtid;
  std::map<int, unsigned long> faccpmtid_to_channelkey_data;

  /// \brief Reset all variables. 
  void ResetVariables();	
 	
  /// \brief ROOT TFile that will be used to store the output from this tool
  TFile* fOutput_tfile = nullptr;

  /// \brief TTree that will be used to store output
  TTree* fPhaseIITrigTree = nullptr;
  TTree* fPhaseIITankClusterTree = nullptr;
  TTree* fPhaseIIMRDClusterTree = nullptr;
 
  std::map<double,std::vector<Hit>>* m_all_clusters = nullptr;  
  std::map<double,std::vector<MCHit>>* m_all_clusters_MC = nullptr;  
  std::map<double,std::vector<unsigned long>>* m_all_clusters_detkeys = nullptr;  
  Geometry *geom = nullptr;

  /// \brief Branch variables
  /// \brief ANNIE event number
  uint32_t fEventNumber;
  uint32_t fRunNumber;
  uint32_t fSubrunNumber;
  uint64_t fEventTimeTank;
  TimeClass* mrd_timestamp=nullptr;
  TimeClass fEventTimeMRD;
  int fRunType;
  uint64_t fStartTime;
  int fNumEntries;
  ULong64_t fEventTimeTank_Tree;
  ULong64_t fEventTimeMRD_Tree;
  ULong64_t fStartTime_Tree;
  int fExtended;	//extended window variable, 0: no extended readout, 1: CC extended readout, 2: Non-CC extended readout
  double fPot;
  int fBeamok;	//1: beam is ok, 0: beam is not ok 
  float BOGUS = 9999.;
  int BOGUS_INT = 9999;
  int BOGUS_INDEX = -1;
  float LOW_FLOAT = -1e30;

  // \brief Event Status flag masks
  int fEventStatusApplied;
  int fEventStatusFlagged;

  // SiPM Hit Info
  int fSiPM1NPulses;
  int fSiPM2NPulses;
  std::vector<double> fSiPMHitQ;
  std::vector<double> fSiPMHitT;
  std::vector<double> fSiPMHitAmplitude;
  std::vector<double> fSiPMNum;
  // Digits (Hits)
  int fNHits = 0;
  std::vector<int> fIsFiltered;
  std::vector<double> fHitX;
  std::vector<double> fHitY;
  std::vector<double> fHitZ;
  std::vector<double> fHitT;
  std::vector<double> fHitQ; 
  std::vector<double> fHitPE; 
  std::vector<int> fHitType;
  std::vector<int> fHitDetID;
  std::vector<int> fHitChankey;
  std::vector<int> fHitChankeyMC;
  
  //Digits
  int fNDigitsPMTs = 0;
  int fNDigitsLAPPDs = 0;
  std::vector<double> fdigitX;
  std::vector<double> fdigitY;
  std::vector<double> fdigitZ;
  std::vector<double> fdigitT;
  

  // MRD hit info 
  int fVetoHit;
  std::vector<double> fMRDHitT;
  std::vector<int> fMRDHitDetID;
  std::vector<int> fMRDHitChankey;
  std::vector<int> fMRDHitChankeyMC;
  std::vector<double> fFMVHitT;
  std::vector<int> fFMVHitDetID;
  std::vector<int> fFMVHitChankey;
  std::vector<int> fFMVHitChankeyMC;
  std::map<unsigned long,vector<Hit>>* TDCData=nullptr;
  std::map<unsigned long,vector<MCHit>>* TDCData_MC=nullptr;

  // ************** MRD Cluster level information ********** //
  int fMRDClusterNumber;
  int fMRDClusterHits;
  double fMRDClusterTime;
  double fMRDClusterTimeSigma;
  // Cluster properties
  std::vector<double> mrddigittimesthisevent;
  std::vector<int> mrddigitpmtsthisevent;
  std::vector<unsigned long> mrddigitchankeysthisevent;
  std::vector<std::vector<int>> MrdTimeClusters;
  
  // ************** Tank Cluster level information ********** //
  std::map<double,double> ClusterMaxPEs;
  std::map<double,Position> ClusterChargePoints;
  std::map<double,double> ClusterChargeBalances;
  int fClusterNumber;
  int fNumClusterTracks;
  int fClusterHits;
  double fClusterCharge;
  double fClusterTime;
  double fClusterPE;
  double fClusterMaxPE;
  double fClusterChargePointX;
  double fClusterChargePointY;
  double fClusterChargePointZ;
  double fClusterChargeBalance;
  std::vector<int> fADCWaveformChankeys; 
  std::vector<int> fADCWaveformSamples;  

  // ************ Muon reconstruction level information ******** //
  std::string MRDTriggertype;
  std::vector<BoostStore>* theMrdTracks;   // the reconstructed tracks
  int numtracksinev;
  std::vector<double> fMRDTrackAngle;
  std::vector<double> fMRDTrackAngleError;
  std::vector<double> fMRDPenetrationDepth;
  std::vector<double> fMRDTrackLength;
  std::vector<double> fMRDEntryPointRadius;
  std::vector<double> fMRDEnergyLoss;
  std::vector<double> fMRDEnergyLossError;
  std::vector<double> fMRDTrackStartX;
  std::vector<double> fMRDTrackStartY;
  std::vector<double> fMRDTrackStartZ;
  std::vector<double> fMRDTrackStopX;
  std::vector<double> fMRDTrackStopY;
  std::vector<double> fMRDTrackStopZ;
  std::vector<bool> fMRDSide;
  std::vector<bool> fMRDStop;
  std::vector<bool> fMRDThrough;
  int fMRDContained;

  // Trigger-level information
  std::map<std::string,bool> fDataStreams;
  int fTriggerword;
  int fTankMRDCoinc;
  int fNoVeto;
  int fHasTank;
  int fHasMRD;

  // ************ MC Truth Information **************** //
  uint64_t fMCEventNum;
  uint16_t fMCTriggerNum;
  int fiMCTriggerNum;
  // True muon
  double fTrueVtxX;
  double fTrueVtxY;
  double fTrueVtxZ;
  double fTrueVtxTime;
  double fTrueDirX;
  double fTrueDirY;
  double fTrueDirZ;
  double fTrueAngle;
  double fTruePhi;
  double fTrueMuonEnergy;
  double fTrueMuonMomentum;
  int fTruePrimaryPdg;
  double fTrueTrackLengthInWater; 
  double fTrueTrackLengthInMRD; 
  std::vector<int> *fTruePrimaryPdgs = nullptr;
  std::vector<double> *fTrueNeutCapVtxX = nullptr;
  std::vector<double> *fTrueNeutCapVtxY = nullptr;
  std::vector<double> *fTrueNeutCapVtxZ = nullptr;
  std::vector<double> *fTrueNeutCapNucleus = nullptr;
  std::vector<double> *fTrueNeutCapTime = nullptr;
  std::vector<double> *fTrueNeutCapGammas = nullptr;
  std::vector<double> *fTrueNeutCapE = nullptr;
  std::vector<double>* fTrueNeutCapGammaE = nullptr;
  int fTrueMultiRing;
  double fTrueDistanceToEdge;

  //Weights
  std::map<std::string, std::vector<double>> fxsec_weights;
  std::map<std::string, std::vector<double>> fflux_weights;
  std::vector<double> fAll;
  std::vector<double> fAxFFCCQEshape;
  std::vector<double> fDecayAngMEC;
  std::vector<double> fNormCCCOH;
  std::vector<double> fNorm_NCCOH;
  std::vector<double> fRPA_CCQE;
  std::vector<double> fRootinoFix;
  std::vector<double> fThetaDelta2NRad;
  std::vector<double> fTheta_Delta2Npi;
  std::vector<double> fTunedCentralValue;
  std::vector<double> fVecFFCCQEshape;
  std::vector<double> fXSecShape_CCMEC;
  std::vector<double> fpiplus;
  std::vector<double> fpiminus;
  std::vector<double> fkplus;
  std::vector<double> fkzero;
  std::vector<double> fkminus;
  std::vector<double> fhorncurrent;
  std::vector<double> fpioninexsec;
  std::vector<double> fpionqexsec;
  std::vector<double> fpiontotxsec;
  std::vector<double> fexpskin;
  std::vector<double> fnucleoninexsec;
  std::vector<double> fnucleonqexsec;
  std::vector<double> fnucleontotxsec;

  //Genie information for event
  double fTrueNeutrinoEnergy;
  double fTrueNeutrinoMomentum_X;
  double fTrueNeutrinoMomentum_Y;
  double fTrueNeutrinoMomentum_Z;
  double fTrueNuIntxVtx_X;
  double fTrueNuIntxVtx_Y;
  double fTrueNuIntxVtx_Z;
  double fTrueNuIntxVtx_T;
  double fTrueFSLVtx_X;
  double fTrueFSLVtx_Y;
  double fTrueFSLVtx_Z;
  double fTrueFSLMomentum_X;
  double fTrueFSLMomentum_Y;
  double fTrueFSLMomentum_Z;
  double fTrueFSLTime;
  double fTrueFSLMass;
  int fTrueFSLPdg;
  double fTrueFSLEnergy;
  double fTrueQ2;
  double fTrueW2;
  double fTrueBJx;
  double fTruey;
  double fTrueq0;
  double fTrueq3;
  int fTrueTarget;
  int fTrueCC;
  int fTrueNC;
  int fTrueQEL;
  int fTrueRES;
  int fTrueDIS;
  int fTrueCOH;
  int fTrueMEC;
  int fTrueNeutrons;
  int fTrueProtons;
  int fTruePi0;
  int fTruePiPlus;
  int fTruePiPlusCher;
  int fTruePiMinus;
  int fTruePiMinusCher;
  int fTrueKPlus;
  int fTrueKPlusCher;
  int fTrueKMinus;
  int fTrueKMinusCher;

  // Pion and kaon counts for event
  int fPi0Count;
  int fPiPlusCount;
  int fPiMinusCount;
  int fK0Count;
  int fKPlusCount;
  int fKMinusCount;

/////////////////////////////////////////////////////////////////
// CC0pi Signal Categories 
/////////////////////////////////////////////////////////////////
  const int kUnknown = 0;

  // Signal events broken down by underlying reaction mode
  const int kSignalCCQE = 1;
  const int kSignalCCMEC = 2;
  const int kSignalCCRES = 3;
  const int kSignalOther = 4;

  // True numu CC event with at least one final-state pion above threshold
  const int kNuMuCCNpi = 5;

  // True numu CC event with zero final-state pions above threshold and
  // zero final-state protons above threshold
  const int kNuMuCC0pi0p = 6;

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  const int kNuMuCCOther = 7;

  // True nue CC event
  const int kNuECC = 8;

  // True neutral current event for any neutrino flavor
  const int kNC = 9;

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  const int kOOFV = 10;

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  const int kOther = 11;
///////////////////////////////////////////////////////////////////////

  // Signal definition requirements
  int fTrueIs_mc;
  int fTruemc_neutrino_is_numu;
  int fTruemc_vertex_in_FV;
  int fTruemc_muon_in_mom_range;
  int fTruemc_lead_p_in_mom_range;
  int fTruemc_no_fs_mesons;
  // Intersection of all of these requirements
  //bool mc_is_signal_ ;

  // Extra flags for looking specifically at final-state pions
  int fTruemc_no_fs_pi0;
  int fTruemc_no_charged_pi_above_threshold;

  int fTruemc_is_cc0pi_signal;

  int fTruemc_num_protons;
  int fTruemc_num_neutrons;
  int fTruemc_num_charged_pions;
  int fTruemc_num_neutral_pions;

  // Water Cherenkov cuts
  int fTrue_category;

  const int ELECTRON_NEUTRINO = 12;
  const int MUON = 13;
  const int MUON_NEUTRINO = 14;
  const int TAU_NEUTRINO = 16;
  const int PROTON = 2212;
  const int NEUTRON = 2112;
  const int PI_ZERO = 111;
  const int PI_PLUS = 211;

  // Checked in GeV 
  const double MUON_P_MIN_MOM_CUT = 600;  // MeV/c
  const double MUON_P_MAX_MOM_CUT = 1200; // MeV/c
  const double MUON_FORWARD_GOING_CUT = 0.8;
  const double CHARGED_PI_MOM_CUT = .160; //water Cherenkov momentum threshold GeV

  const double LEAD_P_MIN_MOM_CUT = 1.100;  // GeV //water Cherenkov momentum threshold of Proton 
  const double LEAD_P_MAX_MOM_CUT = 1.100; // GeV  currently not using a maxium threshold for protons

  // **************** Full reco chain information ************* //
  //  seed vertices
  std::vector<double> fSeedVtxX;
  std::vector<double> fSeedVtxY;
  std::vector<double> fSeedVtxZ;
  std::vector<double> fSeedVtxFOM;
  double fSeedVtxTime;
  
  // Reco vertex
  // Point Position Vertex
  double fPointPosX;
  double fPointPosY;
  double fPointPosZ;
  double fPointPosTime;
  double fPointPosFOM;
  int fPointPosStatus;
  double fPointDirX;
  double fPointDirY;
  double fPointDirZ;
  double fPointDirTime;
  double fPointDirFOM;
  int fPointDirStatus;
  
  // Point Vertex Finder
  double fPointVtxPosX;
  double fPointVtxPosY;
  double fPointVtxPosZ;
  double fPointVtxTime;
  double fPointVtxDirX;
  double fPointVtxDirY;
  double fPointVtxDirZ;
  double fPointVtxFOM;
  int fPointVtxStatus;

  // Simple Reco
  int fSimpleFlag;
  double fSimpleEnergy;
  double fSimpleMomentum;
  double fSimpleMomentumCor;
  double fSimpleVtxX;
  double fSimpleVtxY;
  double fSimpleVtxZ;
  double fSimpleStopVtxX;
  double fSimpleStopVtxY;
  double fSimpleStopVtxZ;
  double fSimpleCosTheta;
  double fSimplePt;
  int fSimpleFV;
  double fSimpleMrdEnergyLoss;
  double fSimpleTrackLengthInMRD;
  double fSimpleTrackLengthInTank;
  double fSimpleMRDStartX;
  double fSimpleMRDStartY;
  double fSimpleMRDStartZ;
  double fSimpleMRDStopX;
  double fSimpleMRDStopY;
  double fSimpleMRDStopZ;
  double fSimpleDistanceToEdge;

  // Ring Counting
  double fRCSRPred;
  double fRCMRPred;

  // Extended Vertex
  double fRecoVtxX;
  double fRecoVtxY;
  double fRecoVtxZ;
  double fRecoVtxTime;
  double fRecoDirX;
  double fRecoDirY;
  double fRecoDirZ;
  double fRecoVtxFOM;
  double fRecoAngle;
  double fRecoPhi;
  int fRecoStatus;
  
  // ************* Difference between MC and Truth *********** //
  double fDeltaVtxX; 
  double fDeltaVtxY;
  double fDeltaVtxZ;
  double fDeltaVtxR;
  double fDeltaVtxT;
  double fDeltaParallel;
  double fDeltaPerpendicular;
  double fDeltaAzimuth;
  double fDeltaZenith;  
  double fDeltaAngle;
  
  // MuonFitter vertex
  int fJulieRecoFlag;
  double fJulieRecoEnergy;
  double fJulieRecoVtxX;
  double fJulieRecoVtxY;
  double fJulieRecoVtxZ;
  double fJulieRecoStopVtxX;
  double fJulieRecoStopVtxY;
  double fJulieRecoStopVtxZ;
  double fJulieRecoCosTheta;
  double fJulieRecoPt;
  int fJulieRecoFV;
  double fJulieRecoMrdEnergyLoss;
  double fJulieRecoTrackLengthInMRD;
  double fJulieRecoMRDStartVtxX;
  double fJulieRecoMRDStartVtxY;
  double fJulieRecoMRDStartVtxZ; 
  double fJulieRecoMRDStopVtxX;
  double fJulieRecoMRDStopVtxY;
  double fJulieRecoMRDStopVtxZ; 
  double fJulieRecoNeutrinoEnergy;
  double fJulieRecoQ2;
  double fRecoMuonVtxX;
  double fRecoMuonVtxY;
  double fRecoMuonVtxZ;
  double fRecoTankTrack;
  double fRecoMuonKE;
  int fNumMrdLayers;

  /// \brief Integer that determines the level of logging to perform
  int verbosity = 0;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;		
  int get_ok;	

  /// \Integer flags that control additional output to the PhaseIITree
  bool TankClusterProcessing = 0;
  bool MRDClusterProcessing = 0;
  bool TriggerProcessing = 1;
  bool TankHitInfo_fill = 0;
  bool MRDHitInfo_fill = 0;
  bool fillCleanEventsOnly = 0; //Only output events not flagged by EventSelector tool
  bool MCTruth_fill = 0; //Output the MC truth information
  bool TankReco_fill = 0;
  bool MRDReco_fill = 0;
  bool Reweight_fill = 0;
  bool SimpleReco_fill = 0;
  bool RingCounting_fill = 0;
  bool RecoDebug_fill = 0; //Outputs results of Reconstruction at each step (best fits, FOMs, etc.)
  bool muonTruthRecoDiff_fill = 0; //Output difference in tmuonruth and reconstructed values
  bool SiPMPulseInfo_fill = 0;
  bool Digit_fill = 0;
  bool MuonFitter_fill = 0; //juju
  bool VertexDistancetoEdge_fill = 0;
  bool CC0PiSignal_fill = 0; 
};


#endif
