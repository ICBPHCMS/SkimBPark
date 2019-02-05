//author: george karathanasis, georgios.karathanasis@cern.ch
//inspired by https://github.com/CMSBParking/cmssw/blob/872669fa1dd7d9bfb8c97f4ffd961a8a6a43761a/RecoEgamma/EgammaElectronProducers/plugins/LowPtGsfElectronSeedProducer.cc (rob/mauro)
#ifndef EPROPAGATION_H
#define EPROPAGATION_H
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "FWCore/Framework/interface/ESHandle.h"

  using namespace reco;


class Epropagation{
  
 public:
  Epropagation();

  Epropagation(const std::shared_ptr<reco::Track> track_ref,const edm::Handle<reco::PFClusterCollection>& clusters,edm::ESHandle<MagneticField> field_ );

  virtual ~Epropagation();

  bool success() { return m_success; }
 
  float deta(){ return m_deta; }
 
 float dphi(){ return m_dphi; }

 float ep() { return m_ep;}

 private:
//    typedef reco::PFClusterCollection PFClusters; 
    bool m_success;
    float m_deta;
    float m_dphi;
    float  m_ep;


};
#endif
