#include "HTauTauTree.h"
#include <iostream>
#include<math.h> 
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "OfflineProducerHelper.h"
#include "../interface/GenFlags.h"

//For trigger selction tools, please refer to 
//https://github.com/LLRCMS/LLRHiggsTauTau/blob/master/NtupleProducer/Utils/OfflineProducerHelper.h

using namespace std;


//return a certain bit in the gen flags word
bool GenFlag (int word, GenFlags::genFlags bitPos)
{
    int bit = static_cast<int> bitPos;
    bool res = word & ( 1 << bit );
    return res;
}

// not supported now, a dedicated branch was used
/*
GenFlags::genFlags GetHZDecayMode (int word)
{
    if (GenFlags (word), GenFlags::HZToMuHad)   return GenFlags::HZToMuHad;
    if (GenFlags (word), GenFlags::HZToEHad)    return GenFlags::HZToEHad;
    if (GenFlags (word), GenFlags::HZToHadHad)  return GenFlags::HZToHadHad;
    if (GenFlags (word), GenFlags::HZToMuMu)    return GenFlags::HZToMuMu;
    if (GenFlags (word), GenFlags::HZToEE)      return GenFlags::HZToEE;
    if (GenFlags (word), GenFlags::HZToEMu)     return GenFlags::HZToEMu;
    if (GenFlags (word), GenFlags::HZToEEPrompt)   return GenFlags::HZToEEPrompt;
    if (GenFlags (word), GenFlags::HZToMuMuPrompt) return GenFlags::HZToMuMuPrompt;
    if (GenFlags (word), GenFlags::HZToOther)      return GenFlags::HZToOther;

    // default, should never happen if input is a good particle (H/Z)
    return GenFlags::HZToOther;
}
*/

void PlotTreeVars()
{
    //TFile* fIn = new TFile ("../test/HTauTauAnalysis.root"); 
    TFile* fIn = new TFile ("/opt/sbg/scratch1/cms/lebihan/chineseSchool/HiggsTauTau/HTauTauAnalysis.root"); 
    //TFile* fIn = new TFile ("/opt/sbg/scratch1/cms/lebihan/chineseSchool/ZTauTau/HTauTauAnalysis.root"); 
    
    TTree* treePtr = (TTree*) fIn->Get("HTauTauTree/HTauTauTree");
    TH1F *evCounter = (TH1F*) fIn->Get("HTauTauTree/Counters");
    
    cout<<"Number of Events analysed "<<evCounter->GetBinContent(1)<<endl;
    cout<<"Number of Events Passing at least one trigger "<<evCounter->GetBinContent(2)<<endl;
    
    //HTauTauTree class defines the branches from the Ntuplizer output
    HTauTauTree* tree = new HTauTauTree (treePtr);
    
    //trigger helper 
    OfflineProducerHelper helper;
    /*
    // Can activate only some branches like this
    tree->GetTree()->SetBranchStatus("*", 0);
    tree->GetTree()->SetBranchStatus("bDiscriminator", 1);
    tree->GetTree()->SetBranchStatus("mothers_px", 1);
    */

    // histograms
    TH1F* bDiscr = new TH1F ("h_bDiscr", "bDiscriminator", 20, 0, 10);
    TH1F* h_mothers_px = new TH1F ("h_mothers_px", "mothers_px", 100, -100, 100);
    TH1F* h_mothers_SVmass = new TH1F ("h_mothers_SVfitMass", "mothers_SVfitMass", 100, -100, 800);
    TH2F *h_dauPtvsSVFit = new TH2F("dauPtvsSVMass","dauPtvsSVMass",100,-100,800,100,0,1000);
    TH1F* h_mothers_Vismass = new TH1F ("h_mothers_VisMass", "mothers_VisMass", 100, -100, 800);
    TH1F* h_mothers_Colinearmass = new TH1F ("h_mothers_ColinearMass", "mothers_ColinearMass", 100, -100, 800);
    
    //long int nentries = tree->GetEntries();
    long int nentries = 2000;
    //long int nentries = 10;
    cout << "TreeEntries: " << nentries << endl;
    for (int iEntry = 0; iEntry < nentries ; iEntry++)
    {
        //cout << "  -------------------- NEW EVENT --------------------- " << endl;
        
        tree->GetEntry(iEntry);

        //Check if a specific trigger is fired
        if(helper.IsTriggerFired(tree->triggerbit,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1")){
            cout<<"trigger HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1 fired in event "<<iEntry<<endl;
        }
        //print all fired triggers
        //cout<<"NEW EVENT"<<endl<<endl;
        //helper.printFiredPaths(tree->triggerbit);
        
        //Distribution of the bDiscriminant
        for (int iJet = 0; iJet < tree->bDiscriminator->size(); iJet++)
        {
            bDiscr->Fill (tree->bDiscriminator->at(iJet));
        }
        //std::cout << "========" << std::endl;
        //Distributions for the pairs and for the daughters
        for (int iMoth = 0; iMoth < tree->mothers_px->size(); iMoth++)
        {
            h_mothers_px->Fill(tree->mothers_px->at(iMoth));
	  
            float svmass=tree->SVfitMass->at(iMoth);
	    if (svmass != -999.)
	    {
            h_mothers_SVmass->Fill(svmass);
	    //std::cout << "svmass " << svmass << std::endl;
            //Distribution of daughters pt as a function of the pair invariant mass
            //indexes of the daughters (in the daughters vector)
            int dau1index = tree->indexDau1->at(iMoth);
            int dau2index = tree->indexDau2->at(iMoth);
            //daughters four momenta
            float dau1px = tree->daughters_px->at(dau1index);
            float dau1py = tree->daughters_py->at(dau1index);
            float dau1pz = tree->daughters_pz->at(dau1index);
            float dau1e = tree->daughters_e->at(dau1index);
            float dau2px = tree->daughters_px->at(dau2index);
            float dau2py = tree->daughters_py->at(dau2index);
            float dau2pz = tree->daughters_pz->at(dau2index);
            float dau2e = tree->daughters_e->at(dau2index);
            TLorentzVector dau1(dau1px, dau1py, dau1pz, dau1e);
            TLorentzVector dau2(dau2px, dau2py, dau2pz, dau2e);
            h_dauPtvsSVFit->Fill(svmass, dau1.Pt());
            h_dauPtvsSVFit->Fill(svmass, dau2.Pt());
			float vismass=sqrt((dau1e+dau2e)*(dau1e+dau2e)-((dau1px+dau2px)*(dau1px+dau2px)+(dau1py+dau2py)*(dau1py+dau2py)+(dau1pz+dau2pz)*(dau1pz+dau2pz)));
			//cout << vismass << endl;
			
			float det=sin(dau1.Theta())*cos(dau1.Phi())*sin(dau2.Theta())*sin(dau2.Phi())-sin(dau2.Theta())*cos(dau2.Phi())*sin(dau1.Theta())*sin(dau1.Phi());
			
			//float metx=
			//float mety=
			
			
			float mu1p=(metx*sin(dau2.Theta())*sin(dau2.Phi())-sin(dau2.Theta())*cos(dau2.Phi())*mety)/det;
			float mu2p=(sin(dau1.Theta())*cos(dau1.Phi())*mety-metx*sin(dau1.Theta())*sin(dau1.Phi()))/det;
			
			
			//asume same direction
			float mu1px=mu1p*sin(dau1.Theta())*cos(dau1.Phi());
			float mu1py=mu1p*sin(dau1.Theta())*sin(dau1.Phi());
			float mu1pz=mu1p*cos(dau1.Theta());
			
			
			float mu2px=mu2p*sin(dau2.Theta())*cos(dau2.Phi());
			float mu2py=mu2p*sin(dau2.Theta())*sin(dau2.Phi());
			float mu2pz=mu2p*cos(dau2.Theta());
			
			
			float colinearmass=sqrt((dau1e+dau2e+mu1p+mu2p)*(dau1e+dau2e+mu1p+mu2p)-((dau1px+dau2px+mu1px+mu2px)*(dau1px+dau2px+mu1px+mu2px)+(dau1py+dau2py+mu1py+mu2py)*(dau1py+dau2py+mu1py+mu2py)+(dau1pz+dau2pz+mu1pz+mu2pz)*(dau1pz+dau2pz+mu1pz+mu2pz)));
			
			
			h_mothers_Vismass->Fill(vismass);
	        h_mothers_Colinearmass->Fill(colinearmass);
	    
	    //float det = sin(dau1.Theta())...;
	
	   }
        }
        
        // reco leptons
        //cout << " === All saved leptons ===" << endl;        
        for (int iDau = 0; iDau < tree->daughters_px->size(); iDau++)
        {
            float px = tree->daughters_px->at(iDau);        
            float py = tree->daughters_py->at(iDau);
            float pt = sqrt (px*px + py*py);
            //cout << iDau << ") id: " << tree->PDGIdDaughters->at(iDau) << " -- pt :" << pt <<  " | gen: " << tree->daughters_genindex->at(iDau) << endl;
        }
        //cout << endl;
        
        // gen level info
        /*cout << " === All saved gen info ===" << endl;
        for (int iGen = 0; iGen < tree->genpart_px->size(); iGen++)
        {
            int ApdgId = abs (tree->genpart_pdg->at(iGen));
            cout << iGen << ")  Id: " << tree->genpart_pdg->at(iGen);
            if (GenFlag (tree->genpart_flags->at(iGen) , GenFlags::fromH ) )   cout << "   , fromH num. " << tree->genpart_HMothInd->at(iGen);
            if (GenFlag (tree->genpart_flags->at(iGen) , GenFlags::fromTop ) ) cout << "   , fromTop num. " << tree->genpart_TopMothInd->at(iGen);
            if (GenFlag (tree->genpart_flags->at(iGen) , GenFlags::fromTau ) ) cout << "   , fromTau num. " << tree->genpart_TauMothInd->at(iGen);
            if (ApdgId == 25 || ApdgId == 23)
                cout << "  | decayMode: " << tree->genpart_HZDecayMode->at(iGen) << endl; // meaning of this values is in GenHelper.h
            cout << endl;
        }*/
        //cout << endl;
        
    }
    
    TCanvas * c1 = new TCanvas;
    bDiscr->Draw();
    
    TCanvas * c2 = new TCanvas;
    h_mothers_px->Draw();
    
    TCanvas * c3 = new TCanvas();
    h_dauPtvsSVFit->Draw("COLZ");
    
    TCanvas *c4 = new TCanvas();
    //h_mothers_SVmass->Draw();
    h_mothers_Vismass->Draw("");
    //h_mothers_SVmass->Draw("same");
    h_mothers_Colinearmass->Draw("same");
    c4->SaveAs("c4.pdf");
   
    
    //c2->Print ("mothers_px.png", "png");
    
}

// if compiling: c++ -lm -o PlotTreeVars PlotTreeVars.C `root-config --glibs --cflags`
int main()
{PlotTreeVars();}
