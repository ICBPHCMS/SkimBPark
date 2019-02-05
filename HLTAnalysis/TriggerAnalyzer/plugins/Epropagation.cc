#include "HLTAnalysis/TriggerAnalyzer/plugins/Epropagation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/Math/interface/LorentzVector.h"

Epropagation::Epropagation():
  m_success(false),m_deta(-99),m_dphi(-99),m_ep(-99){}



Epropagation::Epropagation(const std::shared_ptr<reco::Track> track_ref,const edm::Handle<reco::PFClusterCollection>& clusters,  edm::ESHandle<MagneticField> field_ ){
// Propagate 'electron' to ECAL surface
  float energy = sqrt( pow(0.000511,2.) + track_ref->outerMomentum().Mag2() );
  XYZTLorentzVector mom= XYZTLorentzVector( track_ref->outerMomentum().x(),
					     track_ref->outerMomentum().y(),
					     track_ref->outerMomentum().z(),
					     energy );
  XYZTLorentzVector pos=XYZTLorentzVector( track_ref->outerPosition().x(),
					     track_ref->outerPosition().y(),
					     track_ref->outerPosition().z(),
					     0. );
  // XYZTLorentzVector mom= XYZTLorentzVector( 0,0,0,0 );
  //  XYZTLorentzVector pos=XYZTLorentzVector( 0,0,0,0. );
  math::XYZVector field(field_->inTesla(GlobalPoint(0,0,0)));
  BaseParticlePropagator particle( RawParticle(mom,pos), 0, 0, field.z() );
  particle.setCharge(track_ref->charge());
  particle.propagateToEcalEntrance(false);
  if ( particle.getSuccess() != 0 ) { m_success=false; m_deta=-99; m_dphi=-99; m_ep=-99; }
  else {
  // ECAL entry point for track
  GlobalPoint ecal_pos(particle.vertex().x(), particle.vertex().y(), particle.vertex().z());
  // Preshower limit
  bool below_ps = pow(ecal_pos.z(),2.) > pow(2.50746495928f,2.)*ecal_pos.perp2();

  // Store info for PreId
  struct Info {
    reco::PFClusterRef clu_ref = reco::PFClusterRef();
    float min_dr = 1.e6;
    float dr = 1.e6;
    float deta = 1.e6;
    float dphi = 1.e6;
    math::XYZPoint shower_pos = math::XYZPoint(0.,0.,0.);
  } info;
  
  // Iterate through ECAL clusters 
  for ( unsigned int iclu = 0; iclu < clusters.product()->size(); iclu++ ) {
    reco::PFClusterRef clu_ref(clusters,iclu);

    // Correct ecal_pos for shower depth 
    double shower_depth = reco::PFCluster::getDepthCorrection(clu_ref->correctedEnergy(), below_ps,false);
    GlobalPoint shower_pos = ecal_pos + 
      GlobalVector(particle.momentum().x(),
		   particle.momentum().y(),
		   particle.momentum().z()).unit() * shower_depth;

    // Determine deta, dphi, dr
    float deta = std::abs( clu_ref->positionREP().eta() - shower_pos.eta() );
    float dphi = std::abs( clu_ref->positionREP().phi() - shower_pos.phi() );
    if ( dphi > float(TMath::Pi())) { dphi -= float(TMath::TwoPi()); }
    float dr = std::sqrt( std::pow(dphi,2.f) + std::pow(deta,2.f) );
    
    // Find nearest ECAL cluster
    if ( dr < info.min_dr ) {
      info.min_dr = dr;
      info.clu_ref = clu_ref;
      info.dr = dr;
      info.deta = deta;
      info.dphi = dphi;
      info.shower_pos = shower_pos;
    }
  
  }

  // Populate PreId object
  math::XYZPoint point( ecal_pos.x(),ecal_pos.y(),ecal_pos.z() );
  float ep = info.clu_ref->correctedEnergy() / std::sqrt( track_ref->innerMomentum().mag2() );
  // float ep=0;
  m_deta=info.deta; m_dphi= info.dphi; m_ep=ep; m_success=true;
  }
}


Epropagation::~Epropagation() {}
