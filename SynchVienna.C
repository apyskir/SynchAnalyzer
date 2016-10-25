#define SynchVienna_cxx
#include "SynchVienna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void SynchVienna::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SynchVienna.C
//      Root > SynchVienna t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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
      
      hMuPt->Fill(pt_1);
      hTauPt->Fill(pt_2);
      hMuPhi->Fill(phi_1);
      hTauPhi->Fill(phi_2);
      hMuEta->Fill(eta_1);
      hTauEta->Fill(eta_2);
      hMET->Fill(met);
      hMETphi->Fill(metphi);
      hMud0->Fill(d0_1);
      hTaud0->Fill(d0_2);
      hMudZ->Fill(dZ_1);
      hTaudZ->Fill(dZ_2);
      hMuIso->Fill(iso_1);
      hMuIsoZoom->Fill(iso_1);
      hTauIso->Fill(iso_2);
      hMuMt->Fill(mt_1);
      hTauMt->Fill(mt_2);
      
   }
}
