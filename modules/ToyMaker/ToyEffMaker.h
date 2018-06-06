#ifndef TOY_EFF_MAKER_H
#define TOY_EFF_MAKER_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"


#include "TRandom3.h"
#include "TLorentzVector.h"


#include "vendor/loguru.h"

class ToyEffMaker : public HistoAnalyzer {
protected:
	XmlFunction xfPt, xfEta, xfPhi, xfpPt, xfpM;
	size_t N = 0;


	CutCollection dcuts;

	map< string, XmlHistogram> xh;
	map< string, TH3*> h3;

public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		// loguru::add_file( TString::Format("ToyMaker_%d.log", config.getInt("jobIndex") ).Data(), loguru::Truncate, loguru::Verbosity_MAX);

		LOG_F( INFO, "pt distribution @ %s", config.q( "dist.TF1{name==Pt}" ).c_str() );
		xfPt.set( config, config.q( "dist.TF1{name==Pt}" ) );
		LOG_F( INFO, "Eta distribution @ %s", config.q( "dist.TF1{name==Eta}" ).c_str() );
		xfEta.set( config, config.q( "dist.TF1{name==Eta}" ) );
		LOG_F( INFO, "Phi distribution @ %s", config.q( "dist.TF1{name==Phi}" ).c_str() );
		xfPhi.set( config, config.q( "dist.TF1{name==Phi}" ) );

		LOG_F( INFO, "pairPt distribution @ %s", config.q( "dist.TF1{name==pairPt}" ).c_str() );
		if ( config.exists( config.q( "dist.TF1{name==pairPt}" ) ) )
			xfpPt.set( config, config.q( "dist.TF1{name==pairPt}" ) );
		
		LOG_F( INFO, "pairMass distribution @ %s", config.q( "dist.TF1{name==pairMass}" ).c_str() );
		xfpM.set( config, config.q( "dist.TF1{name==pairMass}" ) );


		N = config.get<size_t>( "N" );

		XmlHistogram txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==mc_5}" ) ); xh["mc_5"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==mc_6}" ) ); xh["mc_6"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==rc_5}" ) ); xh["rc_5"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==rc_6}" ) ); xh["rc_6"] = txh;

		txh.load( config, config.q( nodePath + ".XmlHistogram{name==pos_eff}" ) ); xh["pos_mtd_eff"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==neg_eff}" ) ); xh["neg_mtd_eff"] = txh;

		for ( auto kv : xh ){
			h3[ kv.first ] = (TH3*)xh[ kv.first ].getTH1().get();
		}
		
		assert( h3["rc_5"] != nullptr );
		assert( h3["rc_6"] != nullptr );
		assert( h3["mc_5"] != nullptr );
		assert( h3["mc_6"] != nullptr );

		assert( h3["pos_mtd_eff"] != nullptr );
		assert( h3["neg_mtd_eff"] != nullptr );

		h3[ "rc_eff_5" ] = (TH3*)h3[ "rc_5" ]->Clone("rc_eff_5");
		h3[ "rc_eff_6" ] = (TH3*)h3[ "rc_6" ]->Clone("rc_eff_6");

		h3[ "rc_eff_5" ]->Divide( h3[ "mc_5" ] );
		h3[ "rc_eff_6" ]->Divide( h3[ "mc_6" ] );


		for ( auto kv : h3 ){
			LOG_F( INFO, "H3[%s] = %p", kv.first.c_str(), kv.second );
		}

		gRandom = new TRandom3();
		unsigned long long int seed = get_seed();
		LOG_F( INFO, "seed=%llu", seed );
		gRandom->SetSeed( seed );


		LOG_F( INFO, "Loading cuts on daughter kinematics" );
		dcuts.init( config, "p.daughter" );
		vector<string> cutNames = dcuts.names();
		for ( string n : cutNames ){
			LOG_F( INFO, "DaughterCut [ %s ] = %s", n.c_str(), dcuts[ n ]->toString().c_str()  );
		}

	}

	unsigned long long int get_seed(){
		unsigned long long int random_value = 0; //Declare value to store data into
		size_t size = sizeof(random_value); //Declare size of data
		ifstream urandom("/dev/urandom", ios::in|ios::binary); //Open stream
		if(urandom) //Check if stream is open
		{
			urandom.read(reinterpret_cast<char*>(&random_value), size); //Read from urandom
			if(urandom) {
				return random_value;
			}
			else { //Read failed
				return 0;
			}
			urandom.close(); //close stream
		} else { //Open failed
			std::cerr << "Failed to open /dev/urandom" << std::endl;
		}
		return 0;
	}

	TLorentzVector gen(){
		TLorentzVector lv;
		
		float pt = xfPt.getTF1()->GetRandom();
		float eta = xfEta.getTF1()->GetRandom();
		float phi = xfPhi.getTF1()->GetRandom();
		while ( pt < 0 || fabs( eta ) > 0.5 ){
			pt = xfPt.getTF1()->GetRandom();
		}
		lv.SetPtEtaPhiM( pt, eta, phi, 0.105 );
		return lv;
	}


	double effRC( double pt, double eta, double phi, int charge ){
		TH3 * h = h3[ "rc_eff_5" ];
		if ( -1 == charge )
			h = h3[ "rc_eff_6" ];
		int binx = h->GetXaxis()->FindBin( phi );
		int biny = h->GetYaxis()->FindBin( eta );
		int binz = h->GetZaxis()->FindBin( pt );
		int bin = h->GetBin( binx, biny, binz );
		LOG_F( 1, "(pt, eta, phi) = (%f, %f, %f)", pt, eta, phi );
		LOG_F( 1, "(bx, by, bz) = (%d, %d, %d) = %d", binx, biny, binz, bin );
		LOG_F( 1, "v = %f", h->GetBinContent( bin )  );
		return h->GetBinContent( bin );
	}
	double effRC( TLorentzVector &lv, int charge ){
		return effRC( lv.Pt(), lv.PseudoRapidity(), lv.Phi(), charge );
	} 


	double effMTD( double pt, double eta, double phi, int charge ){
		TH3 * h = h3[ "pos_mtd_eff" ];
		if ( -1 == charge )
			h = h3[ "neg_mtd_eff" ];
		int binx = h->GetXaxis()->FindBin( phi );
		int biny = h->GetYaxis()->FindBin( eta );
		int binz = h->GetZaxis()->FindBin( pt );
		int bin = h->GetBin( binx, biny, binz );
		LOG_F( 1, "(pt, eta, phi) = (%f, %f, %f)", pt, eta, phi );
		LOG_F( 1, "(bx, by, bz) = (%d, %d, %d) = %d", binx, biny, binz, bin );
		LOG_F( 1, "v = %f", h->GetBinContent( bin )  );
		return h->GetBinContent( bin );
	}
	double effMTD( TLorentzVector &lv, int charge ){
		return effMTD( lv.Pt(), lv.PseudoRapidity(), lv.Phi(), charge );
	} 


	virtual void applyBoost( TLorentzVector &_parent_lv, TLorentzVector &_d_lv ){

		float betaX = _parent_lv.Px() / _parent_lv.E();
		float betaY = _parent_lv.Py() / _parent_lv.E();
		float betaZ = _parent_lv.Pz() / _parent_lv.E();

		_d_lv.Boost(betaX,betaY,betaZ);
	}


	/* Two Body Decay in the rest frame of the parent
	 *
	 * ONLY Valid for m1 == m2! The daughter mass must be equal for this simplified form.
	 * The decay is computed and then the daughters are boosted into the frame of the parent
	 */
	virtual void twoBodyDecay( TLorentzVector _parent_lv, TLorentzVector &daughter1, TLorentzVector &daughter2 ){
		// DEBUG( classname(), "lv( P=" << dts(_parent_lv.Px()) << "," << dts(_parent_lv.Py()) << "," << dts(_parent_lv.Pz()) << ", M=" << dts( _parent_lv.M() ) << ")" );

		// ParticleInfo &d1 = products[0];
		// ParticleInfo &d2 = products[1];

		// MUST BE EQUAL
		double M_d1 = 0.105;
		double M_d2 = 0.105;


		double E_d = _parent_lv.M()/2.;
		double p = sqrt(E_d*E_d - M_d1*M_d1);
		double costheta = gRandom->Uniform(-1.,1.);
		double phi = gRandom->Uniform(0,TMath::Pi()*2);

		// make sure that the magnitude of the mom vector is correct!
		// May allow these distributions to be input?
		double pz = p*costheta;
		double px = p*sqrt(1.-costheta*costheta)*cos(phi);
		double py = p*sqrt(1.-costheta*costheta)*sin(phi);

		// TLorentzVector daughter1( px, py, pz, E_d);
		// TLorentzVector daughter2( -px, -py, -pz, E_d );
		// 
		daughter1.SetPxPyPzE( px, py, pz, E_d );
		daughter2.SetPxPyPzE( -px, -py, -pz, E_d );

		applyBoost( _parent_lv, daughter1 );
		applyBoost( _parent_lv, daughter2 );

		// Note this is slightly different than Bingchu/Shuai's code
		// That code gets d2's P correct but mass wrong
		// this method gets the entire 4-vector correct

		// daughter1;
		// daughter2;
	}

	virtual void makeWithDaughterKinematics(){
		LOG_SCOPE_FUNCTION(INFO);
		LOG_F( INFO, "Generating N=%lu Pairs", N*N );
		LOG_F( INFO, "xfPt = %p", xfPt.getTF1().get() );

		TLorentzVector lv1, lv2, lv;
		size_t nPass = 0;
		for ( size_t i = 0; i < N; i++ ){
			lv1 = gen();
			// lv1.SetPtEtaPhiM( 1.2, 0, -1.0, 0.105 );
			double e1 = effRC( lv1, 1 );
			double emtd = effMTD( lv1, 1 );

			book->fill( "rc_eff_pt", lv1.Pt(), e1 );
			book->fill( "rc_eff_eta", lv1.PseudoRapidity(), e1 );
			book->fill( "rc_eff_phi", lv1.Phi(), e1 );

			book->fill( "mtd_pt", lv1.Pt() );
			book->fill( "mtd_eta", lv1.PseudoRapidity() );
			book->fill( "mtd_phi", lv1.Phi() );

			book->fill( "mtd_pt_w", lv1.Pt(), emtd );
			book->fill( "mtd_eta_w", lv1.PseudoRapidity(), emtd );
			book->fill( "mtd_phi_w", lv1.Phi(), emtd );

			for ( size_t j = 0; j < N; j++ ){
				
				lv2 = gen();

				lv = lv1 + lv2;

				jdb::progressBar( j  + i * N, N*N );

				if ( abs(lv.Rapidity()) > 0.5 )
					continue;


				nPass++;

				double eff = 1.0;
				eff *= e1 * effRC( lv2, -1 );
				book->fill( "hmc", lv.M(), lv.Pt() );
				book->fill( "hrc", lv.M(), lv.Pt(), eff );

				double emtd2 = effMTD( lv2, -1 );
				eff *= emtd * emtd2;
				book->fill( "hmtd", lv.M(), lv.Pt(), eff );
				
			} // loop j
		} // loop i

		LOG_F( INFO, "%lu / %lu = %f", nPass, N*N, nPass / (float)(N*N) );
	}

	virtual void makeWithParentKinematics(){
		LOG_SCOPE_FUNCTION(INFO);
		TLorentzVector lv1, lv2, lv;
		size_t nPass = 0;
		size_t nAttempt = 0;

		while ( nPass < N ){
			nAttempt++;
			// Sample the parent's kinematics
			double m   = xfpM.getTF1()->GetRandom();
			double pt  = xfPt.getTF1()->GetRandom();
			double eta = xfEta.getTF1()->GetRandom();
			double phi = xfPhi.getTF1()->GetRandom();
			lv.SetPtEtaPhiM( pt, eta, phi, m );

			// Require parent to be flat in allowed rapidity range
			while ( abs(lv.Rapidity() ) > 0.5 ){
				eta = xfEta.getTF1()->GetRandom();
				lv.SetPtEtaPhiM( pt, eta, phi, m );
			}

			book->fill( "hgen", lv.M(), lv.Pt() );
			// fill parent kinematic histograms
			book->fill( "gen_pt", lv.Pt() );
			book->fill( "gen_eta", lv.PseudoRapidity() );
			book->fill( "gen_rap", lv.Rapidity() );
			book->fill( "gen_phi", lv.Phi() );
			

			// Perform the 2-body decay
			twoBodyDecay( lv, lv1, lv2 );

			// Daughter kinematic requirement
			// just repeate decay if the eta region is out of range, faster than starting over
			while( 	false == dcuts["eta"]->inInclusiveRange( lv1.PseudoRapidity() ) || 
					false == dcuts["eta"]->inInclusiveRange( lv2.PseudoRapidity() ) ){
				twoBodyDecay( lv, lv1, lv2 );
			}

			//  Daughter Kinematic requirements
			if ( 	false == dcuts["pt"]->inInclusiveRange( lv1.Pt() ) || 
					false == dcuts["pt"]->inInclusiveRange( lv2.Pt() ) )
				continue;

			nPass++;
			jdb::progressBar( nPass, N );

			// fill parent kinematic histograms
			book->fill( "pt", lv.Pt() );
			book->fill( "eta", lv.PseudoRapidity() );
			book->fill( "rap", lv.Rapidity() );
			book->fill( "phi", lv.Phi() );

			double e1 = effRC( lv1, 1 );
			double emtd1 = effMTD( lv1, 1 );

			double e2 = effRC( lv2, -1 );
			double emtd2 = effMTD( lv2, -1 );

			book->fill( "mc_pt", lv1.Pt() );
			book->fill( "mc_eta", lv1.PseudoRapidity() );
			book->fill( "mc_y", lv1.Rapidity() );
			book->fill( "mc_phi", lv1.Phi() );

			book->fill( "mc_pt", lv2.Pt() );
			book->fill( "mc_eta", lv2.PseudoRapidity() );
			book->fill( "mc_eta", lv2.Rapidity() );
			book->fill( "mc_phi", lv2.Phi() );

			//  RC weights only
			book->fill( "rc_pt_w"  , lv1.Pt()             , e1 );
			book->fill( "rc_eta_w" , lv1.PseudoRapidity() , e1 );
			book->fill( "rc_phi_w" , lv1.Phi()            , e1 );

			book->fill( "rc_pt_w"  , lv2.Pt()             , e2 );
			book->fill( "rc_eta_w" , lv2.PseudoRapidity() , e2 );
			book->fill( "rc_phi_w" , lv2.Phi()            , e2 );

			// MTD weights only
			book->fill( "mtd_pt_w"  , lv1.Pt()             , emtd1 );
			book->fill( "mtd_eta_w" , lv1.PseudoRapidity() , emtd1 );
			book->fill( "mtd_phi_w" , lv1.Phi()            , emtd1 );

			book->fill( "mtd_pt_w"  , lv2.Pt()             , emtd2 );
			book->fill( "mtd_eta_w" , lv2.PseudoRapidity() , emtd2 );
			book->fill( "mtd_phi_w" , lv2.Phi()            , emtd2 );

			// RC and MTD weights
			book->fill( "rcmtd_pt_w"  , lv1.Pt()             , emtd1 * e1 );
			book->fill( "rcmtd_eta_w" , lv1.PseudoRapidity() , emtd1 * e1 );
			book->fill( "rcmtd_phi_w" , lv1.Phi()            , emtd1 * e1 );

			book->fill( "rcmtd_pt_w"  , lv2.Pt()             , emtd2 * e2 );
			book->fill( "rcmtd_eta_w" , lv2.PseudoRapidity() , emtd2 * e2 );
			book->fill( "rcmtd_phi_w" , lv2.Phi()            , emtd2 * e2 );


			double eff = 1.0;
			eff *= e1 * e2;
			book->fill( "hmc", lv.M(), lv.Pt() );
			book->fill( "hrc", lv.M(), lv.Pt(), eff );

			eff *= emtd1 * emtd2;
			book->fill( "hmtd", lv.M(), lv.Pt(), eff );


		} // while nPass < N

		LOG_F( INFO, "%lu / %lu = %f", nPass, nAttempt, nPass / (float)(nAttempt) );
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION( INFO );

		book->cd();
		book->makeAll( config, nodePath + ".histograms" );

		if ( "pair" != config["dist:use"] )
			makeWithDaughterKinematics();
		else 
			makeWithParentKinematics();

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
	} // make

};



#endif