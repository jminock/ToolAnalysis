/// This Tool takes the Vertex postion and calcalcates the distance to the tank edge

#include "VtxPointDistanceToTank.h"
double DisplacementTRUE_Z =  168.1; //cm
double DisplacementTRUE_Y =  14.46; //cm


VtxPointDistanceToTank::VtxPointDistanceToTank():Tool(){}


bool VtxPointDistanceToTank::Initialise(std::string configfile, DataModel &data){
  Log("===========================================================================================",v_debug,verbosity);
  Log("VtxPointDistanceToTank Tool: Initialise",v_debug,verbosity);
  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("IsData",isData);
  m_variables.Get("HasGenie",hasGenie);
  m_variables.Get("MakeSimple",makeSimple);
  m_variables.Print();
  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  // THis is Called once 
  //auto get_geometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",geom);
  //if(!get_geometry){
  //	Log("VtxPointDistanceToTank Tool: Error retrieving Geometry from ANNIEEvent!",v_error,verbosity); 
  	//return false; 
  //}

  /////////////
  // Fill Varibles
  //////////////

  //if(makeSimple){
  //this->GetSimpleRECO_Vertex();
  //}
  //
  // if(hasGenie){
  //this->GetTRUE_Vertex();
  //}



  return true;
}

/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////

bool VtxPointDistanceToTank::Execute(){
  Log("===========================================================================================",v_debug,verbosity);
  Log("VtxPointDistanceToTank Tool: Executing",v_debug,verbosity);

   this->ResetVariables();
   if(verbosity>1) ANNIEGeometry::Instance()->PrintGeometry();


if(makeSimple){
  this->GetSimpleRECO_Vertex();
  //std::cout<<" GOT SIMPLE RECO vertex (X,Y,Z)= ("<< fSimpleVtxX*100<< ", "<< fSimpleVtxY*100<< ", "<< fSimpleVtxZ*100<< " ) "<< std::endl;
  fSimpleVtx_DistanceToEdge = ANNIEGeometry::Instance()->DistanceToEdge(fSimpleVtxX*100, fSimpleVtxY*100 + 14.46, fSimpleVtxZ*100 - 168.1);

  Log("fSimpleVtx_DistanceToEdge = " + std::to_string(fSimpleVtx_DistanceToEdge), v_debug, verbosity);
  //if(true){std::cout<<" fSimpleVtx_DistanceToEdge = "<< fSimpleVtx_DistanceToEdge<< std::endl;


  m_data->Stores["RecoEvent"]->Set("SimpleVtx_DistanceToEdge",fSimpleVtx_DistanceToEdge);
}

if(hasGenie){
  this->GetTRUE_Vertex();
   //std::cout<<" GOT TRUE vertex (X,Y,Z)= ("<< fTrueVtxX<< ", "<< fTrueVtxY<< ", "<< fTrueVtxZ<< " ) "<< std::endl;
  fTrue_DistanceToEdge = ANNIEGeometry::Instance()->DistanceToEdge(fTrueVtxX, fTrueVtxY + 14.46, fTrueVtxZ - 168.1);
  //if(true){std::cout<<" fTrue_DistanceToEdge = "<< fTrue_DistanceToEdge<< std::endl;}
    Log("fTrue_DistanceToEdge = " + std::to_string(fTrue_DistanceToEdge), v_debug, verbosity);
    m_data->Stores["GenieInfo"]->Set("True_DistanceToEdge",fTrue_DistanceToEdge);
}





  return true;
}

/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////

bool VtxPointDistanceToTank::Finalise(){

  return true;
}

/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
bool VtxPointDistanceToTank::GetSimpleRECO_Vertex(){

  auto* reco_event = m_data->Stores["RecoEvent"];
  Position SimpleRecoVtx_test;
  if (!reco_event) {
    Log("Error: The VtxPointDistanceToTank tool could not find the RecoEvent Store", v_error, verbosity);
  }

  auto get_vtx = m_data->Stores["RecoEvent"]->Get("SimpleRecoVtx",SimpleRecoVtx_test);

    if( get_vtx){
    fSimpleVtxX = SimpleRecoVtx_test.X();
    fSimpleVtxY = SimpleRecoVtx_test.Y();
    fSimpleVtxZ = SimpleRecoVtx_test.Z();
    //std::cout<<" GOT Simpty vertex  (X,Y,Z)= ("<< fSimpleVtxX<< ", "<< fSimpleVtxY<< ", "<< fSimpleVtxZ<< " ) "<< std::endl;
         Log("GOT (Simple) vertex (X,Y,Z)= ( " + std::to_string(fSimpleVtxX) + ", "  +  std::to_string(fSimpleVtxY) + 
     std::to_string(fSimpleVtxZ) + ")", v_debug, verbosity);


    }
    else{
           Log("Failed to Get Simple RecoVtx ", v_debug, verbosity);
        }

return get_vtx; 
}
/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
bool VtxPointDistanceToTank::GetTRUE_Vertex(){

    //RecoVertex* truevtx = 0;
    //auto get_muonMC = m_data->Stores.at("RecoEvent")->Get("TrueVertex",truevtx);

 double TrueNuIntxVtx_X ,TrueNuIntxVtx_Y ,TrueNuIntxVtx_Z;


    bool get_neutrino_vtxx = m_data->Stores["GenieInfo"]->Get("NuIntxVtx_X",TrueNuIntxVtx_X);
    bool get_neutrino_vtxy = m_data->Stores["GenieInfo"]->Get("NuIntxVtx_Y",TrueNuIntxVtx_Y);
    bool get_neutrino_vtxz = m_data->Stores["GenieInfo"]->Get("NuIntxVtx_Z",TrueNuIntxVtx_Z);
   //auto get_muonMC = m_data->Stores.at("RecoEvent")->Get("TrueVertex",truevtx);

  if (!get_neutrino_vtxx ||
           !get_neutrino_vtxy ||
           !get_neutrino_vtxz) {
    Log("Error: The VtxPointDistanceToTank tool could not find the GetTRUEVertex in RecoEvent Store", v_error, verbosity);
  }

  else if( get_neutrino_vtxx &&
           get_neutrino_vtxy &&
           get_neutrino_vtxz){
    // = truevtx->GetPosition().X();
    // = truevtx->GetPosition().Y();
    // = truevtx->GetPosition().Z();
    fTrueVtxX = TrueNuIntxVtx_X;
    fTrueVtxY = TrueNuIntxVtx_Y;
    fTrueVtxZ = TrueNuIntxVtx_Z;

     Log("GOT vertex (X,Y,Z)= ( " + std::to_string(fTrueVtxX) + ", " +
        std::to_string(fTrueVtxY) + ", " +
        std::to_string(fTrueVtxZ) + ") ", v_debug, verbosity);

     //std::cout<<" GOT vertex (X,Y,Z)= ("<< fTrueVtxX<< ", "<< fTrueVtxY<< ", "<< fTrueVtxZ<< " ) "<< std::endl;
    }
   else{
   //std::cout<<"Failed to Get Simple TrueVtx"<< std::endl;
    Log("Failed to Get Simple TrueVtx", v_debug, verbosity);

   }

return get_neutrino_vtxz&&get_neutrino_vtxx&&get_neutrino_vtxy; 
}
/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////

void VtxPointDistanceToTank::ResetVariables() {

fTrueVtxX= -9999;
fTrueVtxY= -9999;
fTrueVtxZ= -9999;
fTrue_DistanceToEdge= -9999;
fRecoVtxX= -9999;
fRecoVtxY= -9999;
fRecoVtxZ= -9999;
fReco_DistanceToEdge= -9999;
fSimpleVtxX= -9999;
fSimpleVtxY= -9999;
fSimpleVtx_DistanceToEdge= -9999;



}

/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
