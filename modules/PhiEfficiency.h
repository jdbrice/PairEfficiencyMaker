#ifndef PHI_EFFICIENCY_H
#define PHI_EFFICIENCY_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"


#include "TRandom3.h"
#include "TLorentzVector.h"

#include "CintFunctionLibrary.h"

#include "vendor/loguru.h"


#include <stdlib.h>

const double PHI_MASS = 1.019455;




class PhiEfficiency : public HistoAnalyzer {
protected:
	XmlFunction xfPtResolution, xfPhi;
	size_t N = 0;

	

	CutCollection dcuts;

	map< string, XmlHistogram> xh;
	map< string, TH3*> h3;

	XmlHistogram xtbw;
	shared_ptr<TH1> htbw = nullptr;
	shared_ptr<TF1> tbw;
	shared_ptr<TF1> massDistribution;
	shared_ptr<TF1> ceres;
	// shared_ptr<TF1> ptResolution;

	vector<TH1*> mtdResEff;
	TH2 * hModvsEta = nullptr;
	TH2 * hBLvsPhi = nullptr;

	// use the one inside CintFunctionLibrary
	/* CERES rapidity parameterization
	 * y - rapidity
	 * sqrt_s = collision com energy 
	 * m = hadron mass  hmm 
	 */
	// double CERES( double y, double sqrt_s, double m ){
		
	// 	double a = y*y / ( 2 * sqrt_s / m );
	// 	double b = 4 * sigmaL( sqrt_s, 0.939 /* GeV/c^2 */ ) * ( 1 - a );
	// 	double c = cosh( 3*y / b );

	// 	return pow( c,-2.);
	// }
	// double sigmaL( double sqrt_s, double nucleon_mass ){
	// 	return sqrt( log( sqrt_s / (2 * nucleon_mass) ) );
	// }


	double rapidityToEta( double _y, double _pT, double m ){
		double mT = sqrt(_pT*_pT+m*m);
		double pZ = mT*TMath::SinH(_y);
		double pTot = sqrt(_pT*_pT+pZ*pZ);

		double eta = 0.5 * log( ( pTot+pZ )/(pTot-pZ));
		
		return eta;
	}




public:

	PhiEfficiency() {}
	~PhiEfficiency() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();
		book->cd();

		// loguru::add_file( TString::Format("JPsi_%d.log", config.get<int>("jobIndex") ).Data(), loguru::Truncate, loguru::Verbosity_MAX);

		LOG_F( INFO, "Phi distribution @ %s", config.q( "dist.TF1{name==Phi}" ).c_str() );
		LOG_F( INFO, "pt Resolution @ %s", config.q( "dist.TF1{name==ptResolution}" ).c_str() );

		xfPhi.set( config, config.q( "dist.TF1{name==Phi}" ) );

		xfPtResolution.set( config, config.q( "dist.TF1{name==ptResolution}" ) );
		assert( xfPhi.getTF1() );
		assert( xfPtResolution.getTF1() );

		if ( config.get<bool>( "p.writeResolution", false ) ){
			xfPtResolution.getTF1()->Write();
		}

		N = config.get<size_t>( "N" );
		LOG_F( INFO, "N=%lu", N );

		XmlHistogram txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==mc_5}" ) ); xh["mc_5"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==mc_6}" ) ); xh["mc_6"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==rc_5}" ) ); xh["rc_5"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==rc_6}" ) ); xh["rc_6"] = txh;

		txh.load( config, config.q( nodePath + ".XmlHistogram{name==pos_eff}" ) ); xh["pos_mtd_eff"] = txh;
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==neg_eff}" ) ); xh["neg_mtd_eff"] = txh;

		for ( auto kv : xh ){
			h3[ kv.first ] = (TH3*)xh[ kv.first ].getTH1().get();
			assert( h3[ kv.first ]  && kv.first.c_str() );
		}

		txh.load( config, config.q( nodePath + ".XmlHistogram{name==hModvsEta}" ) ); 
		hModvsEta = (TH2*)txh.getTH1().get()->Clone("hModvsEta");
		hModvsEta->SetDirectory(0);
		txh.load( config, config.q( nodePath + ".XmlHistogram{name==hBLvsPhi}" ) ); 
		hBLvsPhi = (TH2*)txh.getTH1().get()->Clone("hBLvsPhi");
		hBLvsPhi->SetDirectory(0);

		LOG_F( INFO, "hModvsEta=%p", hModvsEta );
		LOG_F( INFO, "hBLvsPhi=%p", hBLvsPhi );
		assert( hModvsEta && hBLvsPhi );
		assert( hModvsEta->GetXaxis() );
		assert( hBLvsPhi->GetXaxis() );
		
		// assert( h3["rc_5"] != nullptr );
		// assert( h3["rc_6"] != nullptr );
		// assert( h3["mc_5"] != nullptr );
		// assert( h3["mc_6"] != nullptr );

		// assert( h3["pos_mtd_eff"] != nullptr );
		// assert( h3["neg_mtd_eff"] != nullptr );

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


		xtbw.load( config, config.q( nodePath + ".XmlHistogram{name==TBW}" ) );
		htbw = xtbw.getTH1();

		double phiM =  PHI_MASS;
		double phiW = 4.26e-3;

		
		

		massDistribution = shared_ptr<TF1>( new TF1( "phi_mass", BreitWigner, 0, 10, 2 ) );

		// Set the BreitWigner to use the width and mass of this plc
		massDistribution->SetParameter( 0, phiW );
		massDistribution->SetParameter( 1, phiM );

		massDistribution->SetRange( phiM - 100 * phiW, phiM + 100 * phiW ); // TODO: make configurable
		// absurdly high resolution ( < 1 MeV/bin) but it only needs to compute CDF once since these should be true 1D functions
		massDistribution->SetNpx(500);


		ceres = shared_ptr<TF1>( new TF1( "phi_y", "pow(cosh(3*x/4./sqrt(log([0]/2./[1]))/(1-pow(x,2)/2./[0]*[2])),-2.)" ) );
		// p0 = sqrt_s
		// p1 = nucleon mass
		// p2 = hadron mass
		ceres->SetParameters( 200, 0.938, phiM);
		ceres->SetRange( -1, 1 );


		// force distributions to initialize
		LOG_F( INFO, "initialize MASS" );
		double m   = massDistribution->GetRandom();
		LOG_F( INFO, "initialize CERES" );
		double eta = rapidityToEta( ceres->GetRandom(), 1.0, PHI_MASS );
		// should be flat
		LOG_F( INFO, "initialize Phi" );
		double phi = xfPhi.getTF1()->GetRandom();


		bool tbwMakeHisto = config.get<bool>( "p.tbwMakeHisto", true );

		if( tbwMakeHisto ){
			
			LOG_F( INFO, "building TBW function, this is slow" );
			tbw = shared_ptr<TF1>(  TsallisBlastWave( "phi_tbw", phiM, 0, 0.0964, 0, 1.0926, 40, -0.8, 0.8 )  );
			// tbw->SetNpx(5000);
			tbw->SetRange( 0, 15 );

			assert( tbw );
			LOG_F( INFO, "Makeing TBW Histo so we can speed up future runs" );
			TFile *fTBW = new TFile( "Phi_TBW.root", "RECREATE" );
			TH1*h = new TH1F( "TBW", ";p_{T} (GeV/c); d^{2}N/2#pidp_{T}dp_{T}dy", 15000, 0, 15 );
			for ( size_t iphi = 0; iphi < 1000000000; iphi++ ){
				h->Fill( tbw->GetRandom() );
				if ( iphi % 1000000 == 0 ){
					cout << "." << std::flush;
				}
			}
			cout << endl; 

			fTBW->Write();
			exit(0);
		}

		// load the MTD response Eff tables
		TFile *fMtdResEff = new TFile( config.get<TString>( "p.MtdResEff:url", "Run15ResponseEffViaPtTemplate.root" ) );
		assert( fMtdResEff );
		for ( size_t iBL = 0; iBL < 30; iBL++ ){
			for ( size_t iMod = 0; iMod < 5; iMod++){
				TString h1name = TString::Format( "MtdRespEffCosmic_BL%lu_Mod%lu", iBL+1, iMod+1 );
				TH1 * h1MtdEff = (TH1*)fMtdResEff->Get( h1name );
				mtdResEff.push_back( h1MtdEff );
				if ( nullptr != h1MtdEff ){
					LOG_F( INFO, "Loaded MtdResEff BL%lu, Mod%lu", iBL, iMod );
				} else {
					LOG_F( WARNING, "Cannot Load MtdResEff BL%lu, Mod%lu", iBL, iMod );
				}
			}
		}

	} // initialize

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




	virtual void makeWithParentKinematics(){
		LOG_SCOPE_FUNCTION(INFO);
		TLorentzVector lv1, lv2, lv;
		
		size_t nPass = 0;
		size_t nAttempt = 0;

		while ( nPass < N ){
			nAttempt++;
			
			// Sample the parent's kinematics
			double m   = massDistribution->GetRandom();
			double pt  = htbw->GetRandom();
			double eta = rapidityToEta( ceres->GetRandom(), pt, PHI_MASS );
			// should be flat
			double phi = xfPhi.getTF1()->GetRandom();
			
			lv.SetPtEtaPhiM( pt, eta, phi, m );

			// Require parent to be flat in allowed rapidity range
			while ( abs(lv.Rapidity() ) > 0.5 ){
				eta = rapidityToEta( ceres->GetRandom(), pt, PHI_MASS );
				lv.SetPtEtaPhiM( pt, eta, phi, m );
			}

			// LOG_F( INFO, "GEN: m=%0.3f, pT=%0.3f, eta=%0.3f, phi=%0.3f", m, pt, eta, phi );

			book->fill( "hgen", lv.M(), lv.Pt() );
			// fill parent kinematic histograms
			book->fill( "gen_pt", lv.Pt() );
			book->fill( "gen_eta", lv.PseudoRapidity() );
			book->fill( "gen_rap", lv.Rapidity() );
			book->fill( "gen_phi", lv.Phi() );
			

			// Perform the 2-body decay
			twoBodyDecay( lv, lv1, lv2 );

			double rcpt1 = smearPt( lv1.Pt() );
			double rcpt2 = smearPt( lv2.Pt() );

			TLorentzVector rclv, rclv1, rclv2;
			LOG_F( 1, "pt1 Smeared %f->%f", lv1.Pt(), rcpt1 );
			LOG_F( 1, "pt2 Smeared %f->%f", lv2.Pt(), rcpt2 );
			rclv1.SetPtEtaPhiM( rcpt1, lv1.Eta(), lv1.Phi(), 0.105 );
			rclv2.SetPtEtaPhiM( rcpt2, lv2.Eta(), lv2.Phi(), 0.105 );
			rclv = rclv1 + rclv2;

			// LOG_F( INFO, "lv1(pt=%0.3f, eta=%0.3f, phi=%0.3f)", lv1.Pt(), lv1.Eta(), lv1.Phi() );
			// LOG_F( INFO, "lv2(pt=%0.3f, eta=%0.3f, phi=%0.3f)", lv2.Pt(), lv2.Eta(), lv2.Phi() );

			nPass++;
			jdb::progressBar( nPass, N );
			// continue;
			// Daughter kinematic requirement
			// just repeate decay if the eta region is out of range, faster than starting over
			
			//  Daughter Kinematic requirements
			if ( 	false == dcuts["pt"]->inInclusiveRange( rclv1.Pt() ) || 
					false == dcuts["pt"]->inInclusiveRange( rclv2.Pt() ) ){
				continue;
			}

			if ( 	false == dcuts["eta"]->inInclusiveRange( rclv1.Eta() ) || 
					false == dcuts["eta"]->inInclusiveRange( rclv2.Eta() ) ){
				continue;
			}


			book->fill( "hmc", rclv.M(), rclv.Pt() );
			// fill parent kinematic histograms
			book->fill( "mc_pt", rclv.Pt() );
			book->fill( "mc_eta", rclv.PseudoRapidity() );
			book->fill( "mc_rap", rclv.Rapidity() );
			book->fill( "mc_phi", rclv.Phi() );

			// intentionally use the MC info for efficiency lookup since tables are made using MC kinematics
			double e1 = effRC( lv1, 1 );
			double emtd1 = effMTD( lv1, 1 );
			double emtdRes1 = effMTDRes( lv1.Pt(), lv1.Eta(), lv1.Phi() );

			double e2 = effRC( lv2, -1 );
			double emtd2 = effMTD( lv2, -1 );
			double emtdRes2 = effMTDRes( lv2.Pt(), lv2.Eta(), lv2.Phi() );

			book->fill( "dmc_pt", rclv1.Pt() );
			book->fill( "dmc_eta", rclv1.PseudoRapidity() );
			book->fill( "dmc_y", rclv1.Rapidity() );
			book->fill( "dmc_phi", rclv1.Phi() );

			book->fill( "dmc_pt", rclv2.Pt() );
			book->fill( "dmc_eta", rclv2.PseudoRapidity() );
			book->fill( "dmc_eta", rclv2.Rapidity() );
			book->fill( "dmc_phi", rclv2.Phi() );

			//  RC weights only
			book->fill( "drc_pt_w"  , rclv1.Pt()             , e1 );
			book->fill( "drc_eta_w" , rclv1.PseudoRapidity() , e1 );
			book->fill( "drc_phi_w" , rclv1.Phi()            , e1 );

			book->fill( "drc_pt_w"  , rclv2.Pt()             , e2 );
			book->fill( "drc_eta_w" , rclv2.PseudoRapidity() , e2 );
			book->fill( "drc_phi_w" , rclv2.Phi()            , e2 );

			// MTD weights only
			book->fill( "dmtd_pt_w"  , rclv1.Pt()             , emtd1 );
			book->fill( "dmtd_eta_w" , rclv1.PseudoRapidity() , emtd1 );
			book->fill( "dmtd_phi_w" , rclv1.Phi()            , emtd1 );

			book->fill( "dmtd_pt_w"  , rclv2.Pt()             , emtd2 );
			book->fill( "dmtd_eta_w" , rclv2.PseudoRapidity() , emtd2 );
			book->fill( "dmtd_phi_w" , rclv2.Phi()            , emtd2 );

			// RC and MTD weights
			book->fill( "drcmtd_pt_w"  , rclv1.Pt()             , emtd1 * e1 );
			book->fill( "drcmtd_eta_w" , rclv1.PseudoRapidity() , emtd1 * e1 );
			book->fill( "drcmtd_phi_w" , rclv1.Phi()            , emtd1 * e1 );

			book->fill( "drcmtd_pt_w"  , rclv2.Pt()             , emtd2 * e2 );
			book->fill( "drcmtd_eta_w" , rclv2.PseudoRapidity() , emtd2 * e2 );
			book->fill( "drcmtd_phi_w" , rclv2.Phi()            , emtd2 * e2 );


			double eff = 1.0;
			eff *= e1 * e2;
			book->fill( "hrc", rclv.M(), rclv.Pt(), eff );

			book->fill( "rc_pt", rclv.Pt(), eff );
			book->fill( "rc_eta", rclv.PseudoRapidity(), eff );
			book->fill( "rc_rap", rclv.Rapidity(), eff );
			book->fill( "rc_phi", rclv.Phi(), eff );

			eff *= emtd1 * emtd2;
			book->fill( "hmtd", rclv.M(), rclv.Pt(), eff );

			book->fill( "mtd_pt", rclv.Pt(), eff );
			book->fill( "mtd_eta", rclv.PseudoRapidity(), eff );
			book->fill( "mtd_rap", rclv.Rapidity(), eff );
			book->fill( "mtd_phi", rclv.Phi(), eff );


			eff *= emtdRes1 * emtdRes2;
			book->fill( "mtdr_pt", rclv.Pt(), eff );
			book->fill( "mtdr_eta", rclv.PseudoRapidity(), eff );
			book->fill( "mtdr_rap", rclv.Rapidity(), eff );
			book->fill( "mtdr_phi", rclv.Phi(), eff );


		} // while nPass < N

		LOG_F( INFO, "%lu / %lu = %f", nPass, nAttempt, nPass / (float)(nAttempt) );
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION( INFO );

		book->cd();
		book->makeAll( config, nodePath + ".histograms" );

		makeWithParentKinematics();

		if ( 0 == config.get<int>( "jobIndex", -2 ) || -1 == config.get<int>( "jobIndex", -2 ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
	} // make


	double effRC( double pt, double eta, double phi, int charge ){
		TH3 * h = h3[ "rc_eff_5" ];
		if ( -1 == charge )
			h = h3[ "rc_eff_6" ];
		int binx = h->GetXaxis()->FindBin( phi );
		int biny = h->GetYaxis()->FindBin( eta );
		int binz = h->GetZaxis()->FindBin( pt );
		int bin = h->GetBin( binx, biny, binz );
		

		if ( binx >= h->GetXaxis()->GetNbins() || binx == 0 )
			return 0;
		if ( biny >= h->GetYaxis()->GetNbins() || biny == 0 )
			return 0;
		if ( binz >= h->GetZaxis()->GetNbins() || binz == 0 )
			return 0;


		if ( eta > 1.0 ){
			LOG_F( INFO, "(pt, eta, phi) = (%f, %f, %f)", pt, eta, phi );
			LOG_F( INFO, "(bx, by, bz) = (%d, %d, %d) = %d nBins( %d, %d, %d )", binx, biny, binz, bin, h->GetXaxis()->GetNbins(), h->GetYaxis()->GetNbins(), h->GetZaxis()->GetNbins() );
			LOG_F( INFO, "v = %f", h->GetBinContent( bin )  );
		}

		

		return h->GetBinContent( bin );
	}
	double effRC( TLorentzVector &lv, int charge ){
		return effRC( lv.Pt(), lv.PseudoRapidity(), lv.Phi(), charge );
	} 


	double smearPt(double pt){
		TF1 * momResolution = xfPtResolution.getTF1().get();
		// =============================== MOMENTUM SMEARING ===============================
		double ptRes = momResolution->Eval( pt );
		LOG_F( 1, "ptRes = %0.3f @ pT = %0.3f", ptRes, pt );
		double rndCrystalBall = gRandom->Gaus( 0, 1.0 );
		return pt * (1 + rndCrystalBall * ptRes  );
	}


	double effMTDRes( double pt, double eta, double phi ) {

		// project the hModvsEta and sample for a Mod
		int iEta1 = hModvsEta->GetXaxis()->FindBin( eta - 0.05 );
		int iEta2 = hModvsEta->GetXaxis()->FindBin( eta + 0.05 );
		TH1 * hModAtEta = hModvsEta->ProjectionY( "tempMod", iEta1, iEta2 );
		int iMod = floor( hModAtEta->GetRandom() );
		if ( iMod < 0 || iMod > 4 ){
			LOG_F( ERROR, "mod=%d @ eta = %f", iMod, eta );
			return 0;
		}
		

		// project the hBLvsPhi and sample for a BL
		int iPhi1 = hBLvsPhi->GetXaxis()->FindBin( phi - 0.05 );
		int iPhi2 = hBLvsPhi->GetXaxis()->FindBin( phi + 0.05 );
		TH1 * hBLAtPhi = hBLvsPhi->ProjectionY( "tempBL", iPhi1, iPhi2 );
		int iBL = floor( hBLAtPhi->GetRandom() );
		if ( iBL < 0 || iBL > 29 ){
			LOG_F( INFO, "bl=%d @ phi = %f", iBL, phi );
			return 0;
		}

		// generate a random BL / Mod
		int i = 5 * iBL + iMod;
		if ( i >= mtdResEff.size() ){
			assert(0 && "overflow");
		}
		TH1 * h1 = mtdResEff[ i ];
		int bin = h1->GetXaxis()->FindBin( pt );
		if ( bin == 0 || bin >= h1->GetXaxis()->GetNbins() )
			return 0;
		double r = h1->GetBinContent( bin );
		if ( r < 0 || r > 1){
			LOG_F(ERROR, "r=%f @ pt=%f, iBL=%d, iMod=%d", r, pt, iBL, iMod);
			r = 0;
		} else {
			// LOG_F( INFO,  "r=%f @ pt=%f, iBL=%d, iMod=%d", r, pt, iBL, iMod);
		}

		if ( pt < 1.3 ){
			r *= 500;
			if ( r > 1.0)
				r = 1.0;
		} else {
			r *= 1.05;
		}

		return r;
	}

	double effMTD( double pt, double eta, double phi, int charge ){
		TH3 * h = h3[ "pos_mtd_eff" ];
		if ( -1 == charge )
			h = h3[ "neg_mtd_eff" ];
		int binx = h->GetXaxis()->FindBin( phi );
		int biny = h->GetYaxis()->FindBin( eta );
		int binz = h->GetZaxis()->FindBin( pt );

		if ( binx >= h->GetXaxis()->GetNbins() || binx == 0 )
			return 0;
		if ( biny >= h->GetYaxis()->GetNbins() || biny == 0 )
			return 0;
		if ( binz >= h->GetZaxis()->GetNbins() || binz == 0 )
			return 0;
		
		int bin = h->GetBin( binx, biny, binz );
		LOG_F( 1, "(pt, eta, phi) = (%f, %f, %f)", pt, eta, phi );
		LOG_F( 1, "(bx, by, bz) = (%d, %d, %d) = %d nBins( %d, %d, %d )", binx, biny, binz, bin, h->GetXaxis()->GetNbins(), h->GetYaxis()->GetNbins(), h->GetZaxis()->GetNbins() );
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

};



#endif