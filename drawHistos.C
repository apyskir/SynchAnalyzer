#include "Synch.C"
#include "SynchCERN.C"
#include "SynchLIP.C"
#include "TCanvas.h"

void draw1(TH1F *hWAW, TH1F *hVienna, bool abs){

  TCanvas *c1 = new TCanvas(hWAW->GetName(),"",700,700);

  std::string hName(hWAW->GetName());
  std::string refName(hVienna->GetName());
  std::string hTitle = hName.substr(1, hName.find("_")-1);
  std::string checkName = hName.substr(hName.find("_")+1,hName.length() - hName.find("_"));
  refName = refName.substr(refName.find("_")+1,refName.length() - refName.find("_"));
  c1->SetTitle(hTitle.c_str());
  
  Float_t nWAW = hWAW->Integral(1,hWAW->GetNbinsX()), nVienna = hVienna->Integral(1,hVienna->GetNbinsX());
  
  hWAW->SetStats(kFALSE);
  hVienna->SetLineColor(2);
  
  TPad *pad1;
  TPad *pad2;
  TPad *pad3;
  
  //if abs=true the absolute difference between histos is plotted, if abs=false a relative difference is plotted
  if(!abs){
	  c1->Divide(3);
	  pad1 = (TPad*)c1->GetPad(1);
	  pad2 = (TPad*)c1->GetPad(2);
	  pad3 = (TPad*)c1->GetPad(3);
	  
	  pad1->SetPad(0.01,0.50,0.99,0.99);
	  pad2->SetPad(0.01,0.25,0.99,0.50);
	  pad3->SetPad(0.01,0.03,0.99,0.25);
	  pad2->SetFillStyle(4000);
	  pad3->SetFillStyle(4000);
	  }
  else {
  	  c1->Divide(2);
	  pad1 = (TPad*)c1->GetPad(1);
	  pad2 = (TPad*)c1->GetPad(2);
	  
	  pad1->SetPad(0.01,0.30,0.99,0.99);
	  pad2->SetPad(0.01,0.03,0.99,0.29);
	  pad2->SetFillStyle(4000);
  }
  
  pad1->cd();
  
  float percentage = 100*(nVienna - nWAW)/nVienna;
  hWAW->SetTitle(hTitle.c_str());
  
  std::string tmpName(c1->GetTitle());
  tmpName += (" diff " + checkName + " - "+refName);
  
  TH1F *htmp = new TH1F(tmpName.c_str(),tmpName.c_str(),hWAW->GetNbinsX(),hWAW->GetXaxis()->GetXmin(),hWAW->GetXaxis()->GetXmax());
  for(Int_t i = 1; i < htmp->GetNbinsX()+1; i ++){
  	if(!abs) {if(hVienna->GetBinContent(i)!=0) htmp->SetBinContent(i, hWAW->GetBinContent(i)/hVienna->GetBinContent(i)-1); else htmp->SetBinContent(i, 0);}
  		else htmp->SetBinContent(i, hWAW->GetBinContent(i) - hVienna->GetBinContent(i));
  }
  
  tmpName += (" diff " + refName + " - "+checkName);
  TH1F *htmp2 = new TH1F(tmpName.c_str(),tmpName.c_str(),hWAW->GetNbinsX(),hWAW->GetXaxis()->GetXmin(),hWAW->GetXaxis()->GetXmax());
  for(Int_t i = 1; i < htmp2->GetNbinsX()+1; i ++){
  	if(hWAW->GetBinContent(i)!=0) htmp2->SetBinContent(i, hVienna->GetBinContent(i)/hWAW->GetBinContent(i)-1); else htmp2->SetBinContent(i, 0);
  }
  tmpName += (" diff " + checkName + " - "+refName);
  //////////////////////////////////////////////////
  if(hName=="hMuM_WAW") {hWAW->SetMinimum(0.5); pad1->SetLogy();}
  
  hWAW->SetYTitle("Events");
  hWAW->SetXTitle(hTitle.c_str());
  if(hTitle.find("Zoom")!=std::string::npos) hWAW->SetXTitle(hTitle.substr(0,hTitle.find("Zoom")-1).c_str());
  hWAW->GetYaxis()->SetTitleOffset(0.6);
  hWAW->GetYaxis()->SetTitleSize(0.08);
  if(hName.find("Veto")==std::string::npos) {
  	hWAW->Scale(1/hWAW->Integral(1,hWAW->GetNbinsX()+1));
  	hVienna->Scale(1/hVienna->Integral(1,hVienna->GetNbinsX()+1));
  	}
  
  float max = hWAW->GetMaximum();
  if(max < hVienna->GetMaximum()) max = hVienna->GetMaximum();
  float min = hWAW->GetMinimum();
  if(min > hVienna->GetMinimum()) min = hVienna->GetMinimum();
  if(min != 0 && min - 0.05*(max-min) > 0) hWAW->SetMinimum(min - 0.05*(max-min)); else hWAW->SetMinimum(0);
  hWAW->SetMaximum(max + 0.05*(max-min));
  hWAW->Draw();
  hVienna->Draw("same");
  
  TLegend l(0.45,0.62,0.85,0.72,NULL,"brNDC");
  l.SetTextSize(0.04);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  if(percentage>0) l.AddEntry(hWAW,TString::Format("%s: %.0f = %s - %.3f %%",checkName.c_str(),nWAW,refName.c_str(),percentage).Data());
  	else if(percentage<0) l.AddEntry(hWAW,TString::Format("%s: %.0f = %s + %.3f %%",checkName.c_str(),nWAW,refName.c_str(),-percentage).Data());
  	else l.AddEntry(hWAW,TString::Format("%s: %.0f",checkName.c_str(),nWAW).Data());
  l.AddEntry(hVienna,TString::Format("%s: %.0f",refName.c_str(), nVienna).Data());
  l.Draw("same");
  
  pad2->cd();
  
  htmp->SetStats(kFALSE);
  if(!abs) htmp->SetYTitle(TString::Format("#frac{N_{%s} - N_{%s}}{N_{%s}}",checkName.c_str(),refName.c_str(),refName.c_str()).Data());
  	else htmp->SetYTitle(TString::Format("N_{%s} - N_{%s}",checkName.c_str(),refName.c_str()).Data());
  htmp->SetXTitle(hTitle.c_str());
  htmp->GetYaxis()->SetTitleOffset(0.6);
  htmp->GetXaxis()->SetLabelSize(0.09);
  htmp->GetYaxis()->SetLabelSize(0.09);
  htmp->GetYaxis()->SetTitleSize(0.08);
  htmp->SetTitleSize(0.09);
  htmp->Draw();
  
  if(!abs){
	  pad3->cd();
	  htmp2->SetStats(kFALSE);
	  htmp2->SetYTitle(TString::Format("#frac{N_{%s} - N_{%s}}{N_{%s}}",refName.c_str(),checkName.c_str(),checkName.c_str()).Data());
	  htmp2->SetXTitle(hTitle.c_str());
	  htmp2->GetYaxis()->SetTitleOffset(0.6);
	  htmp2->GetXaxis()->SetLabelSize(0.09);
	  htmp2->GetYaxis()->SetLabelSize(0.09);
	  htmp2->GetYaxis()->SetTitleSize(0.08);
	  htmp2->SetTitleSize(0.09);
	  htmp2->Draw();
	  }
  
  c1->Print(TString::Format("fig_png/%s_%s_%s.png",hTitle.c_str(),checkName.c_str(),refName.c_str()).Data());
  c1->Close();

}

void drawHistos(){
	
	std::string sample1 = "WAW", sample2 = "DESYVBF";
	//those variables decide if you want to draw pictures and if you want to investigate the sync trees by reading single events
	bool draw = 0, investigate = 1;
	
	//second argument of Synch constructor should not contain "_"
	
	//!!!!!!!!!!To use CERN tree (steggeman) one has to use SynchCERN class instead of Synch (because of Double_t variables instead of Float_t)!!!!!!!!
	
	//In CERN ntuple there's also no iso_2 (tau iso)
	
	TFile *f1, *f2;
	TTree *t1, *t2;
	Synch *s1, *s2;
	SynchCERN *s1c, *s2c;
	SynchLIP *s1l, *s2l;
	
	if(sample1 == "WAW"){
		f1 = new TFile("/dmj/hep/apyskir/scratch/HTauTau/RootAnalysis/HTTAnalysis/RootAnalysis_SynchNTupleMT.root");
		t1 = (TTree*)f1->Get("Summary/tree");
		s1 = new Synch(t1, sample1);
	}
	
	if(sample1 == "Vienna"){
		f1 = new TFile("/scratch/apyskir/Data/Synch/SYNCFILE_SUSYGluGluToHToTauTau_M-160_mt_spring16_reHLT_v4.root");
		t1 = (TTree*)f1->Get("TauCheck");
		s1 = new Synch(t1, sample1);
	}
	
	if(sample1 == "ViennaSM"){
		f1 = new TFile("/scratch/apyskir/Data/Synch/SYNCFILE_SUSYGluGluToHToTauTau_M-160_mt_spring16_v4_SM.root");
		t1 = (TTree*)f1->Get("TauCheck");
		s1 = new Synch(t1, sample1);
	}
	
	if(sample1 == "KIT"){
		f1 = new TFile("/scratch/apyskir/Data/Synch/SUSYGluGluToHToTauTauM160_mt_RunIISpring16MiniAODv2_13TeV_MINIAOD.root");
		t1 = (TTree*)f1->Get("mt/ntuple");
		s1 = new Synch(t1, sample1);
	}
	
	if(sample1 == "CERN"){
		f1 = new TFile("/scratch/apyskir/Data/Synch/mt_sm_v2.root");
		t1 = (TTree*)f1->Get("tree");
		s1c = new SynchCERN(t1, sample1);
	}
	
	if(sample1 == "CERNVBF"){
		f1 = new TFile("/scratch/apyskir/Data/Synch/mt_sm_vbf_v1.root");
		t1 = (TTree*)f1->Get("tree");
		s1c = new SynchCERN(t1, sample1);
		}
	
	if(sample2 == "WAWNEW"){
		f2 = new TFile("/dmj/hep/apyskir/scratch/HTauTau/RootAnalysis/HTTAnalysis/RootAnalysis_SynchNTupleAfter.root");
		t2 = (TTree*)f2->Get("Summary/tree");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "DESYVBF"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/Htt_mt_VBF_SM_sync_v6.root");
		t2 = (TTree*)f2->Get("TauCheck");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "Vienna"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/SYNCFILE_SUSYGluGluToHToTauTau_M-160_mt_spring16_reHLT_v4.root");
		t2 = (TTree*)f2->Get("TauCheck");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "ViennaVBF"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/SYNCFILE_VBFHToTauTau_M125_mt_spring16_v2_SM.root");
		t2 = (TTree*)f2->Get("TauCheck");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "ViennaSM"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/SYNCFILE_SUSYGluGluToHToTauTau_M-160_mt_spring16_v4_SM.root");
		t2 = (TTree*)f2->Get("TauCheck");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "KIT"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/SUSYGluGluToHToTauTauM160_mt_RunIISpring16MiniAODv2_13TeV_MINIAOD.root");
		t2 = (TTree*)f2->Get("mt/ntuple");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "KITVBF"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/Htt_mt_VBFHToTauTauM125_v2.root");
		t2 = (TTree*)f2->Get("ntuple");
		s2 = new Synch(t2, sample2);
	}
	
	if(sample2 == "LIPVBF"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/VBFHToTauTauM125.root");
		t2 = (TTree*)f2->Get("h2taus/sync/mt");
		s2l = new SynchLIP(t2, sample2);
	}
	
	if(sample2 == "LIP"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/SUSYGluGluToHToTauTauM160.root");
		t2 = (TTree*)f2->Get("h2taus/mt");
		s2l = new SynchLIP(t2, sample2);
	}
	
	if(sample2 == "CERN"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/mt_sm_v2.root");
		t2 = (TTree*)f2->Get("tree");
		s2c = new SynchCERN(t2, sample2);
	}
	
	if(sample2 == "CERNVBF"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/mt_sm_vbf_v1.root");
		t2 = (TTree*)f2->Get("tree");
		s2c = new SynchCERN(t2, sample2);
	}
	
	if(sample2 == "UWVBF"){
		f2 = new TFile("/scratch/apyskir/Data/Synch/htt_UW_mt_vbf_laura.root");
		t2 = (TTree*)f2->Get("muTauEventTree/eventTree");
		s2 = new Synch(t2, sample2);
	}
	
	if(draw){
		if(sample1.find("CERN")==std::string::npos && sample1.find("LIP")==std::string::npos) s1->Loop();
		if(sample1.find("CERN")!=std::string::npos) s1c->Loop();
		if(sample1.find("LIP")!=std::string::npos) s1l->Loop();
		if(sample2.find("CERN")==std::string::npos && sample2.find("LIP")==std::string::npos) s2->Loop();
		if(sample2.find("LIP")!=std::string::npos) s2l->Loop();
		if(sample2.find("CERN")!=std::string::npos) s2c->Loop();
	
		if(sample2.find("CERN")==std::string::npos && sample2.find("LIP")==std::string::npos){
			draw1(s1->hNPV, s2->hNPV, 0);
			draw1(s1->hNPU, s2->hNPU, 0);
			//muon
			draw1(s1->hMuPt, s2->hMuPt, 1);
			draw1(s1->hMuPhi, s2->hMuPhi, 1);
			draw1(s1->hMuMt, s2->hMuMt, 1);
			draw1(s1->hMuEta, s2->hMuEta, 1);
			draw1(s1->hMud0, s2->hMud0, 1);
			draw1(s1->hMudZ, s2->hMudZ, 1);
			draw1(s1->hMuIso, s2->hMuIso, 1);
			draw1(s1->hMuIsoZoom, s2->hMuIsoZoom, 1);
			draw1(s1->hMuM, s2->hMuM, 0);
			draw1(s1->hMuGenMatch, s2->hMuGenMatch, 1);
			//tau
			draw1(s1->hTauPt, s2->hTauPt, 0);
			draw1(s1->hTauPhi, s2->hTauPhi, 0);
			draw1(s1->hTauMt, s2->hTauMt, 0);
			draw1(s1->hTauEta, s2->hTauEta, 0);
			draw1(s1->hTaud0, s2->hTaud0, 0);
			draw1(s1->hTaudZ, s2->hTaudZ, 0);
			draw1(s1->hTauIso, s2->hTauIso, 0);
			draw1(s1->hTauGenMatch, s2->hTauGenMatch, 1);
			//jets
			draw1(s1->hJet1Pt, s2->hJet1Pt, 0);
			draw1(s1->hJet1Phi, s2->hJet1Phi, 0);
			draw1(s1->hJet1Eta, s2->hJet1Eta, 0);
			draw1(s1->hJet2Pt, s2->hJet2Pt, 0);
			draw1(s1->hJet2Phi, s2->hJet2Phi, 0);
			draw1(s1->hJet2Eta, s2->hJet2Eta, 0);
			draw1(s1->hBJet1Pt, s2->hBJet1Pt, 0);
			draw1(s1->hBJet1Phi, s2->hBJet1Phi, 0);
			draw1(s1->hBJet1Eta, s2->hBJet1Eta, 0);
			draw1(s1->hBJet2Pt, s2->hBJet2Pt, 0);
			draw1(s1->hBJet2Phi, s2->hBJet2Phi, 0);
			draw1(s1->hBJet2Eta, s2->hBJet2Eta, 0);
			//MET
			draw1(s1->hMET, s2->hMET, 0);
			draw1(s1->hMETphi, s2->hMETphi, 0);
			draw1(s1->hMVAMET, s2->hMVAMET, 0);
			draw1(s1->hMVAMETphi, s2->hMVAMETphi, 0);
			draw1(s1->hMVACov00, s2->hMVACov00, 1);
			draw1(s1->hMVACov01, s2->hMVACov01, 1);
			draw1(s1->hMVACov10, s2->hMVACov10, 1);
			draw1(s1->hMVACov11, s2->hMVACov11, 1);
			draw1(s1->hPZetaVis, s2->hPZetaVis, 0);
			draw1(s1->hPZetaMiss, s2->hPZetaMiss, 0);
			draw1(s1->hPFPZetaMiss, s2->hPFPZetaMiss, 0);
			//di-tau
			draw1(s1->hPtTT, s2->hPtTT, 0);
			draw1(s1->hMtTot, s2->hMtTot, 0);
			draw1(s1->hMVis, s2->hMVis, 0);
			draw1(s1->hMSV, s2->hMSV, 0);
			//vbf system
			draw1(s1->hMJJ, s2->hMJJ, 0);
			draw1(s1->hJdEta, s2->hJdEta, 0);
			draw1(s1->hJdPhi, s2->hJdPhi, 0);
			//extra lepton vetos
			draw1(s1->hDiLeptonVeto, s2->hDiLeptonVeto, 0);
			draw1(s1->hExtraElectronVeto, s2->hExtraElectronVeto, 0);
			draw1(s1->hExtraMuonVeto, s2->hExtraMuonVeto, 0);
		}
		
		if(sample1.find("CERN")==std::string::npos && sample1.find("LIP")==std::string::npos && sample2.find("LIP")!=std::string::npos){
			draw1(s1->hNPV, s2l->hNPV, 0);
			draw1(s1->hNPU, s2l->hNPU, 0);
			//muon
			draw1(s1->hMuPt, s2l->hMuPt, 0);
			draw1(s1->hMuPhi, s2l->hMuPhi, 0);
			draw1(s1->hMuMt, s2l->hMuMt, 0);
			draw1(s1->hMuEta, s2l->hMuEta, 0);
			draw1(s1->hMud0, s2l->hMud0, 0);
			draw1(s1->hMudZ, s2l->hMudZ, 0);
			draw1(s1->hMuIso, s2l->hMuIso, 0);
			draw1(s1->hMuIsoZoom, s2l->hMuIsoZoom, 0);
			draw1(s1->hMuM, s2l->hMuM, 0);
			draw1(s1->hMuGenMatch, s2l->hMuGenMatch, 0);
			//tau
			draw1(s1->hTauPt, s2l->hTauPt, 0);
			draw1(s1->hTauPhi, s2l->hTauPhi, 0);
			draw1(s1->hTauMt, s2l->hTauMt, 0);
			draw1(s1->hTauEta, s2l->hTauEta, 0);
			draw1(s1->hTaud0, s2l->hTaud0, 0);
			draw1(s1->hTaudZ, s2l->hTaudZ, 0);
			draw1(s1->hTauIso, s2l->hTauIso, 0);
			draw1(s1->hTauGenMatch, s2l->hTauGenMatch, 0);
			//jets
			draw1(s1->hJet1Pt, s2l->hJet1Pt, 0);
			draw1(s1->hJet1Phi, s2l->hJet1Phi, 0);
			draw1(s1->hJet1Eta, s2l->hJet1Eta, 0);
			draw1(s1->hJet2Pt, s2l->hJet2Pt, 0);
			draw1(s1->hJet2Phi, s2l->hJet2Phi, 0);
			draw1(s1->hJet2Eta, s2l->hJet2Eta, 0);
			//MET
			draw1(s1->hMET, s2l->hMET, 0);
			draw1(s1->hMETphi, s2l->hMETphi, 0);
			draw1(s1->hMVAMET, s2l->hMVAMET, 0);
			draw1(s1->hMVAMETphi, s2l->hMVAMETphi, 0);
			draw1(s1->hMVACov00, s2l->hMVACov00, 0);
			draw1(s1->hMVACov01, s2l->hMVACov01, 0);
			draw1(s1->hMVACov10, s2l->hMVACov10, 0);
			draw1(s1->hMVACov11, s2l->hMVACov11, 0);
			draw1(s1->hPZetaVis, s2l->hPZetaVis, 0);
			draw1(s1->hPZetaMiss, s2l->hPZetaMiss, 0);
			draw1(s1->hPFPZetaMiss, s2l->hPFPZetaMiss, 0);
			//di-tau
			draw1(s1->hPtTT, s2l->hPtTT, 0);
			draw1(s1->hMtTot, s2l->hMtTot, 0);
			draw1(s1->hMVis, s2l->hMVis, 0);
			draw1(s1->hMSV, s2l->hMSV, 0);
			//vbf system
			draw1(s1->hMJJ, s2l->hMJJ, 0);
			draw1(s1->hJdEta, s2l->hJdEta, 0);
			draw1(s1->hJdPhi, s2l->hJdPhi, 0);
			//extra lepton vetos
			draw1(s1->hDiLeptonVeto, s2l->hDiLeptonVeto, 0);
			draw1(s1->hExtraElectronVeto, s2l->hExtraElectronVeto, 0);
			draw1(s1->hExtraMuonVeto, s2l->hExtraMuonVeto, 0);
			}
	
		if(sample1.find("CERN")!=std::string::npos || sample2.find("CERN")!=std::string::npos){
			if(sample1.find("CERN")!=std::string::npos){
				if(sample2.find("LIP")==std::string::npos){
					draw1(s1c->hMuPt, s2->hMuPt, 0);
					draw1(s1c->hTauPt, s2->hTauPt, 0);
					draw1(s1c->hMuPhi, s2->hMuPhi, 0);
					draw1(s1c->hTauPhi, s2->hTauPhi, 0);
					draw1(s1c->hMuMt, s2->hMuMt, 0);
					draw1(s1c->hTauMt, s2->hTauMt, 0);
					draw1(s1c->hMuEta, s2->hMuEta, 0);
					draw1(s1c->hTauEta, s2->hTauEta, 0);
					draw1(s1c->hMud0, s2->hMud0, 0);
					draw1(s1c->hTaud0, s2->hTaud0, 0);
					draw1(s1c->hMudZ, s2->hMudZ, 0);
					draw1(s1c->hTaudZ, s2->hTaudZ, 0);
					draw1(s1c->hMuIso, s2->hMuIso, 0);
					draw1(s1c->hMuIsoZoom, s2->hMuIsoZoom, 0);
					draw1(s1c->hMET, s2->hMET, 0);
					draw1(s1c->hMETphi, s2->hMETphi, 0);
					} else {
					draw1(s1c->hNPV, s2l->hNPV, 0);
					draw1(s1c->hNPU, s2l->hNPU, 0);
					//muon
					draw1(s1c->hMuPt, s2l->hMuPt, 0);
					draw1(s1c->hMuPhi, s2l->hMuPhi, 0);
					draw1(s1c->hMuMt, s2l->hMuMt, 0);
					draw1(s1c->hMuEta, s2l->hMuEta, 0);
					draw1(s1c->hMud0, s2l->hMud0, 0);
					draw1(s1c->hMudZ, s2l->hMudZ, 0);
					draw1(s1c->hMuIso, s2l->hMuIso, 0);
					draw1(s1c->hMuIsoZoom, s2l->hMuIsoZoom, 0);
					draw1(s1c->hMuM, s2l->hMuM, 0);
					draw1(s1c->hMuGenMatch, s2l->hMuGenMatch, 0);
					//tau
					draw1(s1c->hTauPt, s2l->hTauPt, 0);
					draw1(s1c->hTauPhi, s2l->hTauPhi, 0);
					draw1(s1c->hTauMt, s2l->hTauMt, 0);
					draw1(s1c->hTauEta, s2l->hTauEta, 0);
					draw1(s1c->hTaud0, s2l->hTaud0, 0);
					draw1(s1c->hTaudZ, s2l->hTaudZ, 0);
					draw1(s1c->hTauIso, s2l->hTauIso, 0);
					draw1(s1c->hTauGenMatch, s2l->hTauGenMatch, 0);
					//jets
					draw1(s1c->hJet1Pt, s2l->hJet1Pt, 0);
					draw1(s1c->hJet1Phi, s2l->hJet1Phi, 0);
					draw1(s1c->hJet1Eta, s2l->hJet1Eta, 0);
					draw1(s1c->hJet2Pt, s2l->hJet2Pt, 0);
					draw1(s1c->hJet2Phi, s2l->hJet2Phi, 0);
					draw1(s1c->hJet2Eta, s2l->hJet2Eta, 0);
					//MET
					draw1(s1c->hMET, s2l->hMET, 0);
					draw1(s1c->hMETphi, s2l->hMETphi, 0);
					draw1(s1c->hMVAMET, s2l->hMVAMET, 0);
					draw1(s1c->hMVAMETphi, s2l->hMVAMETphi, 0);
					draw1(s1c->hMVACov00, s2l->hMVACov00, 0);
					draw1(s1c->hMVACov01, s2l->hMVACov01, 0);
					draw1(s1c->hMVACov10, s2l->hMVACov10, 0);
					draw1(s1c->hMVACov11, s2l->hMVACov11, 0);
					draw1(s1c->hPZetaVis, s2l->hPZetaVis, 0);
					draw1(s1c->hPZetaMiss, s2l->hPZetaMiss, 0);
					draw1(s1c->hPFPZetaMiss, s2l->hPFPZetaMiss, 0);
					//di-tau
					draw1(s1c->hPtTT, s2l->hPtTT, 0);
					draw1(s1c->hMtTot, s2l->hMtTot, 0);
					draw1(s1c->hMVis, s2l->hMVis, 0);
					draw1(s1c->hMSV, s2l->hMSV, 0);
					//vbf system
					draw1(s1c->hMJJ, s2l->hMJJ, 0);
					draw1(s1c->hJdEta, s2l->hJdEta, 0);
					draw1(s1c->hJdPhi, s2l->hJdPhi, 0);
					//extra lepton vetos
					draw1(s1c->hDiLeptonVeto, s2l->hDiLeptonVeto, 0);
					draw1(s1c->hExtraElectronVeto, s2l->hExtraElectronVeto, 0);
					draw1(s1c->hExtraMuonVeto, s2l->hExtraMuonVeto, 0);
					}
			} else {
				draw1(s1->hNPV, s2c->hNPV, 1);
				draw1(s1->hNPU, s2c->hNPU, 1);
				//muon
				draw1(s1->hMuPt, s2c->hMuPt, 1);
				draw1(s1->hMuPhi, s2c->hMuPhi, 1);
				draw1(s1->hMuMt, s2c->hMuMt, 1);
				draw1(s1->hMuPFMt, s2c->hMuPFMt, 1);
				draw1(s1->hMuEta, s2c->hMuEta, 1);
				draw1(s1->hMud0, s2c->hMud0, 1);
				draw1(s1->hMudZ, s2c->hMudZ, 1);
				draw1(s1->hMuIso, s2c->hMuIso, 1);
				draw1(s1->hMuIsoZoom, s2c->hMuIsoZoom, 0);
				draw1(s1->hMuM, s2c->hMuM, 0);
				draw1(s1->hMuGenMatch, s2c->hMuGenMatch, 1);
				//tau
				draw1(s1->hTauPt, s2c->hTauPt, 1);
				draw1(s1->hTauPhi, s2c->hTauPhi, 1);
				draw1(s1->hTauMt, s2c->hTauMt, 1);
				draw1(s1->hTauPFMt, s2c->hTauPFMt, 1);
				draw1(s1->hTauEta, s2c->hTauEta, 1);
				draw1(s1->hTaud0, s2c->hTaud0, 1);
				draw1(s1->hTaudZ, s2c->hTaudZ, 1);
				draw1(s1->hTauIso, s2c->hTauIso, 1);
				draw1(s1->hTauGenMatch, s2c->hTauGenMatch, 1);
				//jets
				draw1(s1->hJet1Pt, s2c->hJet1Pt, 1);
				draw1(s1->hJet1Phi, s2c->hJet1Phi, 1);
				draw1(s1->hJet1Eta, s2c->hJet1Eta, 1);
				draw1(s1->hJet2Pt, s2c->hJet2Pt, 1);
				draw1(s1->hJet2Phi, s2c->hJet2Phi, 1);
				draw1(s1->hJet2Eta, s2c->hJet2Eta, 1);
				draw1(s1->hBJet1Pt, s2c->hBJet1Pt, 0);
				draw1(s1->hBJet1Phi, s2c->hBJet1Phi, 0);
				draw1(s1->hBJet1Eta, s2c->hBJet1Eta, 0);
				draw1(s1->hBJet2Pt, s2c->hBJet2Pt, 0);
				draw1(s1->hBJet2Phi, s2c->hBJet2Phi, 0);
				draw1(s1->hBJet2Eta, s2c->hBJet2Eta, 0);
				//MET
				draw1(s1->hMET, s2c->hMET, 1);
				draw1(s1->hMETphi, s2c->hMETphi, 1);
				//draw1(s1->hMVAMET, s2c->hMET, 0);
				//draw1(s1->hMVAMETphi, s2c->hMETphi, 0);
				draw1(s1->hMVAMET, s2c->hMVAMET, 1);
				draw1(s1->hMVAMETphi, s2c->hMVAMETphi, 1);
				draw1(s1->hMVACov00, s2c->hMVACov00, 1);
				draw1(s1->hMVACov01, s2c->hMVACov01, 1);
				draw1(s1->hMVACov10, s2c->hMVACov10, 1);
				draw1(s1->hMVACov11, s2c->hMVACov11, 1);
				draw1(s1->hPZetaVis, s2c->hPZetaVis, 1);
				draw1(s1->hPZetaMiss, s2c->hPZetaMiss, 1);
				draw1(s1->hPFPZetaMiss, s2c->hPFPZetaMiss, 1);
				//di-tau
				draw1(s1->hPtTT, s2c->hPtTT, 1);
				draw1(s1->hMtTot, s2c->hMtTot, 1);
				draw1(s1->hMVis, s2c->hMVis, 1);
				draw1(s1->hMSV, s2c->hMSV, 1);
				//vbf system
				draw1(s1->hMJJ, s2c->hMJJ, 1);
				draw1(s1->hJdEta, s2c->hJdEta, 1);
				draw1(s1->hJdPhi, s2c->hJdPhi, 1);
				//extra lepton vetos
				draw1(s1->hDiLeptonVeto, s2c->hDiLeptonVeto, 1);
				draw1(s1->hExtraElectronVeto, s2c->hExtraElectronVeto, 1);
				draw1(s1->hExtraMuonVeto, s2c->hExtraMuonVeto, 1);
			}
		}
	}
	
	Int_t nevt = 0, nevt2 = 0, nevt3 = 0;
	Bool_t is2 = 0;
	if(investigate){
	
		for(Int_t i = 0; i<10000; i++){//s2->fChain->GetEntriesFast(); i++){
			s2->b_evt->GetEntry(i);
			s2->b_met->GetEntry(i);
			is2 = false;
			for(Int_t j = 0; j<s1->fChain->GetEntriesFast(); j++){
				s1->b_evt->GetEntry(j);
				if(s1->evt == s2->evt){
					s1->b_mvamet->GetEntry(j);
					s1->b_njetspt20->GetEntry(j);
					s2->b_njetspt20->GetEntry(i);
					//std::cout<<s1->njetspt20<<" "<<s2->njetspt20<<std::endl;
					//if(s1->njetspt20 == 0 && s2->njetspt20 == 0) {nevt++;}
					if(abs(s1->mvamet - s2->met)>0.0001) {//std::cout<<s1->evt<<": WAW: "<<s1->met<<", CERN: "<<s2c->met<<std::endl; 
						//nevt2++;
						std::cout<<s1->njetspt20<<" "<<s2->njetspt20<<std::endl;
						//if(s1->njetspt20 == 0 && s2->njetspt20 == 0) {nevt3++;}
						//std::cout<<s1->evt<<std::endl;
						//std::cout<<s1->mvamet<<" "<<s2->met<<std::endl;
						}
					}
				}
			//if(!is2) nevt++;
			}
		std::cout<<nevt<<std::endl<<nevt2<<std::endl<<nevt3<<std::endl;
		/*
		for(Int_t i = 0; i<s1->fChain->GetEntriesFast(); i++){
			s1->GetEntry(i);
			if(s1->evt==346755) std::cout<<s1->met<<", "<<s1->metphi<<std::endl;
			}
		for(Int_t i = 0; i<s2c->fChain->GetEntriesFast(); i++){
			s2c->GetEntry(i);
			if(s2c->evt==346755) std::cout<<s2c->met<<", "<<s2c->metphi<<std::endl;
			}*/
	}
	
	std::cout<<"Dziala\n";
	return;
}
