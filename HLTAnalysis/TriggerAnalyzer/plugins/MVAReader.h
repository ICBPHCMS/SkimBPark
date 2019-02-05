//Original author: george karathanais georgios.karathanasis@cern.ch

#ifndef MVAREADER_H
#define MVAREADER_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include <memory>


  class MVAReader {

  public:

    //   HeavyObjectCache( const edm::ParameterSet& );

    MVAReader();

    MVAReader(std::string weight );

    std::unique_ptr<const GBRForest> gbr;    

    virtual ~MVAReader();

    void SetWeight(std::string weight);


   private:
    float trk_pt;
    float trk_eta;
    float trk_phi;
    float trk_p;
    float trk_charge;
    float trk_nhits;
    float trk_high_purity;
    float trk_inp;
    float trk_outp;
    float trk_chi2red;
    float preid_trk_ecal_Deta;
    float preid_trk_ecal_Dphi;
    float preid_e_over_p;
    /*float preid_gsf_dpt;
    float preid_trk_gsf_chiratio;
    float preid_gsf_chi2red;*/

  };


#endif // LowPtGsfElectronSeedHeavyObjectCache_h
