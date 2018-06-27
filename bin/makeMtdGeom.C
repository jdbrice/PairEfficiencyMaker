



void makeMtdGeom(){
	TFile * f = new TFile( "/Users/danielbrandenburg/bnl/local/work/dimuonAna/data/PairDst/pair_dst_pp_DNN_N6.root" );
	TTree * PairDst = (TTree*)f->Get("PairDst");

	TFile * fout = new TFile( "MtdGeomMap.root", "RECREATE" );
	PairDst->Draw( "d1_mModule : d1_mEta>>hModvsEta(200, -1, 1, 5, 0, 5)", "", "colz" );
	PairDst->Draw( "d1_mBackleg : d1_mPhi>>hBLvsPhi(600, -3.2, 3.2, 30, 0, 30)", "", "colz" );

	fout->Write();
}