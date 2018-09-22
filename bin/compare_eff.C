




void compare_eff(){

	TFile *f0 = new TFile( "all-phi.root" );
	TFile *f1 = new TFile( "out/eff-table-phi-1p3--1.root" );


	double signalBins[] = { 0, 2.2, 2.5, 3.0, 3.5, 4.0, 6.0, 8.0 };
	// TH1 * mtdr0 = ((TH1*)f0->Get("mtdr_pt"))->Rebin( 6, "mtdr_pt_rb", signalBins );
	// TH1 * mc0   = ((TH1*)f0->Get("mc_pt"))->Rebin( 6, "mc_pt_rb", signalBins );
	// TH1 * gen0   = ((TH1*)f0->Get("gen_pt"))->Rebin( 6, "gen_pt_rb", signalBins );

	// TH1 * mtdr1 = ((TH1*)f1->Get("mtdr_pt"))->Rebin( 6, "mtdr_pt_rb", signalBins );
	// TH1 * mc1   = ((TH1*)f1->Get("mc_pt"))->Rebin( 6, "mc_pt_rb", signalBins );
	// TH1 * gen1   = ((TH1*)f1->Get("gen_pt"))->Rebin( 6, "gen_pt_rb", signalBins );


	TH1 * mtdr0 = (TH1*)f0->Get("mtdr_pt");
	TH1 * mc0   = (TH1*)f0->Get("mc_pt");
	TH1 * gen0   = (TH1*)f0->Get("gen_pt");

	TH1 * mtdr1 = (TH1*)f1->Get("mtdr_pt");
	TH1 * mc1   = (TH1*)f1->Get("mc_pt");
	TH1 * gen1   = (TH1*)f1->Get("gen_pt");

	



	mtdr0->Divide( gen0 );
	mtdr1->Divide( gen1 );


	mtdr0->SetLineColor(kRed);
	mtdr1->SetLineColor(kBlue);

	mtdr0->Draw();
	mtdr1->Draw("same");

	TH1 * hrel = (TH1*)mtdr1->Clone( "hrel" );
	hrel->Divide( mtdr0 );

	TCanvas * can0 = new TCanvas( "can", "" );
	can0->cd();
	hrel->Draw();

}