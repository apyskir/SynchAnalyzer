#define SynchLIP_cxx
#include "SynchLIP.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void SynchLIP::Loop()
{
//   In a ROOT session, you can do:
//      root> .L SynchLIP.C
//      root> SynchLIP t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      //pileup
      hNPV->Fill(npv);
      hNPU->Fill(npu);
      //muon
      hMuPt->Fill(pt_1);
      hMuPhi->Fill(phi_1);
      hMuEta->Fill(eta_1);
      hMud0->Fill(d0_1);
      hMudZ->Fill(dZ_1);
      hMuIso->Fill(iso_1);
      hMuIsoZoom->Fill(iso_1);
      hMuMt->Fill(mt_1);
      hMuM->Fill(m_1);
      hMuGenMatch->Fill(gen_match_1);
      //tau
      hTauPt->Fill(pt_2);
      hTauPhi->Fill(phi_2);
      hTauEta->Fill(eta_2);
      hTaud0->Fill(d0_2);
      hTaudZ->Fill(dZ_2);
      hTauIso->Fill(iso_2);
      hTauMt->Fill(mt_2);
      hTauGenMatch->Fill(gen_match_2);
      //jets
      hJet1Pt->Fill(jpt_1);
      hJet1Phi->Fill(jphi_1);
      hJet1Eta->Fill(jeta_1);
      //hJet2Pt->Fill(jpt_2);
      //hJet2Phi->Fill(jphi_2);
      //hJet2Eta->Fill(jeta_2);
      //MET
      hMET->Fill(met);
      hMETphi->Fill(metphi);
      hMVAMET->Fill(mvamet);
      hMVAMETphi->Fill(mvametphi);
      hMVACov00->Fill(mvacov00);
      hMVACov01->Fill(mvacov01);
      hMVACov10->Fill(mvacov10);
      hMVACov11->Fill(mvacov11);
      hPZetaVis->Fill(pzetavis);
      hPZetaMiss->Fill(pzetamiss);
      hPFPZetaMiss->Fill(pfpzetamiss);
      //di-tau
      hPtTT->Fill(pt_tt);
      //hMtTot->Fill(mt_tot);
      hMVis->Fill(m_vis);
      hMSV->Fill(m_sv);
      //vbf system
      //hMJJ->Fill(mjj);
      //hJdEta->Fill(jdeta);
      //hJdPhi->Fill(jdphi);
      //extra lepton vetos
      if(dimuon_veto < 0) hDiLeptonVeto->Fill(0); else hDiLeptonVeto->Fill(dimuon_veto);
      if(extraelec_veto < 0) hExtraElectronVeto->Fill(0); else hExtraElectronVeto->Fill(extraelec_veto);
      if(extramuon_veto < 0) hExtraMuonVeto->Fill(0); else hExtraMuonVeto->Fill(extramuon_veto);
      
   }
}
