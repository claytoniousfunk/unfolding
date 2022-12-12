///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  diffyX
//  Version : 1
//  Author: Clayton Bennett
//  Date : 16 May 2022
//  Notes:
//  	- Calculate the differential jet cross section in pp data
//
//  Version Notes:
// 
//      
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// #include "/home/clayton/Analysis/RooUnfold/src/RooUnfold.h"
// #include "/home/clayton/Analysis/RooUnfold/src/RooFitUnfold.h"
#include "/home/clayton/Analysis/RooUnfold/src/RooUnfoldResponse.h"
#include "/home/clayton/Analysis/RooUnfold/src/RooUnfoldBayes.h"



void diffyX_v1(){

    gSystem->Load("/home/clayton/Analysis/RooUnfold/libRooUnfold");

    TFile *f1, *f2;

    f1 = TFile::Open("~/Analysis/code/skimming/pp_scan/rootFiles/pp_1p0_L3Mu5_16May22.root");

    f2 = TFile::Open("~/Analysis/code/skimming/PYTHIA_scan/rootFiles/PYTHIA_response_noHLT_pthat15_excludeNeutrinos_16May22.root"); // excluded neutrinos file
    //f2 = TFile::Open("~/Analysis/code/skimming/PYTHIA_scan/rootFiles/PYTHIA_response_noHLT_pthat15_20May22.root");

    TH1D *Njet;

    f1->GetObject("h_jetpt",Njet);

    TH1D *vz;

    f1->GetObject("h_vz",vz);

    TH2D *h_response;

    f2->GetObject("h_matchedRecoJetPt_genJetPt",h_response);

    TH1D *h_meas, *h_truth;

    h_meas = (TH1D*) h_response->ProjectionX();
    
    h_truth = (TH1D*) h_response->ProjectionY();

    
    //RooUnfoldResponse response(100,0,100);
    RooUnfoldResponse response(h_meas,h_truth,h_response,"response","my response",0);


    RooUnfoldBayes unfold(&response, h_meas, 1);

    TH1D* h_reco= (TH1D*) unfold.Hunfold();
    

    // rebin
    const int Nedge = 6;
    double ptAxis[Nedge] = {50,60,80,120,200,500};


    TH1D *Njet_R = (TH1D*) Njet->Rebin(Nedge-1,"Njet_R",ptAxis);

    // Divide out bin width, Lint, and Nevent

    

    for(int i = 0; i < Njet_R->GetSize(); i++){
        double NjetBin = Njet_R->GetBinContent(i);
        double NjetBin_err = Njet_R->GetBinError(i);
        double dEta = 3.2;
        double dPt = Njet_R->GetBinWidth(i);
        double Lint = 103.9e9;
        double Nevt = 1.;
        if(dPt!=0.){
            Njet_R->SetBinContent(i,NjetBin/(dPt*dEta*Lint*Nevt));
            Njet_R->SetBinError(i,NjetBin_err/(dPt*dEta*Lint*Nevt));
        }
    }

    TCanvas *c1 = new TCanvas("c1","differential cross section",800,800);
    c1->cd();
    TPad *p1 = new TPad("p1","p1",0,0,1,1);
    p1->Draw();
    p1->cd();
    p1->SetLeftMargin(0.2);
    p1->SetBottomMargin(0.2);
    p1->SetLogy();
    Njet_R->SetTitle("");
    Njet_R->SetStats(0);
    Njet_R->GetXaxis()->SetTitle("#font[52]{p}_{T}^{recoJet} [GeV]");
    Njet_R->GetYaxis()->SetTitle("d^{2}#sigma^{pp}_{jet} / d#font[52]{p}_{T} d#eta [mb/GeV]");

    Njet_R->Draw();


    TCanvas *c2 = new TCanvas("c2","unfolded spectra",800,800);
    TPad *p2 = new TPad("p2","p2",0,0.4,1,1);
    TPad *p2_1 = new TPad("p2_1","p2_1",0,0,1,0.4);
    p2->Draw();
    p2->cd();
    p2->SetLeftMargin(0.2);
    p2->SetBottomMargin(0.0);
    p2->SetLogy();
    h_truth->SetLineColor(kRed);
    h_meas->SetLineColor(kGreen);
    h_reco->SetLineColor(kBlue);

    h_truth->SetMarkerColor(kRed);
    h_meas->SetMarkerColor(kGreen);
    h_reco->SetMarkerColor(kBlue);

    h_truth->SetMarkerStyle(8);
    h_meas->SetMarkerStyle(8);
    h_reco->SetMarkerStyle(8);

    h_truth->GetXaxis()->SetLabelSize(0.0);
    h_truth->GetXaxis()->SetTitleSize(0.0);

    h_truth->SetStats(0);
    h_truth->GetYaxis()->SetTitle("Entries");
    h_truth->SetTitle("");
    h_truth->Draw();
    h_meas->Draw("same");
    h_reco->Draw("same");

    TLegend *leg = new TLegend(0.6,0.3,0.85,0.55);
    leg->SetBorderSize(0.0);
    leg->AddEntry(h_truth,"truth");
    leg->AddEntry(h_meas,"measured");
    leg->AddEntry(h_reco,"unfolded");
    leg->Draw();

    TLatex *la = new TLatex();
    la->SetTextFont(42);
    la->SetTextSize(0.034);
    la->DrawLatexNDC(0.6,0.85,"PYTHIA #sqrt{#font[52]{s}} = 5 TeV");
    la->DrawLatexNDC(0.6,0.78,"Bayesian unfolding via RooUnfold");

    c2->cd();
    p2_1->Draw();
    p2_1->cd();
    p2_1->SetLeftMargin(0.2);
    p2_1->SetBottomMargin(0.25);
    p2_1->SetTopMargin(0.0);

    TH1D *r_reco = (TH1D*) h_reco->Clone("r_reco");
    r_reco->Divide(h_reco,h_truth,1,1,"B");
    r_reco->GetYaxis()->SetTitle("unfold / truth");
    r_reco->GetXaxis()->SetTitle("#font[52]{p}_{T}^{jet} [GeV]");
    r_reco->GetXaxis()->SetLabelSize(0.06);
    r_reco->GetXaxis()->SetTitleSize(0.07);
    r_reco->GetYaxis()->SetLabelSize(0.06);
    r_reco->GetYaxis()->SetTitleSize(0.07);
    r_reco->SetTitle("");
    r_reco->SetStats(0);
    r_reco->GetYaxis()->SetRangeUser(0.98,1.02);

    r_reco->SetLineColor(kBlack);
    r_reco->SetMarkerColor(kBlack);
    r_reco->SetMarkerStyle(24);
    r_reco->SetMarkerSize(0.8);
    r_reco->Draw();

    TLine *l1 = new TLine(0,1,500,1);
    l1->SetLineStyle(7);
    l1->Draw();

    TCanvas *c3 = new TCanvas("c3","c3",800,800);
    c3->cd();
    TPad *p3 = new TPad("p3","p3",0,0,1,1);
    p3->SetLeftMargin(0.2);
    p3->SetBottomMargin(0.2);
    p3->SetRightMargin(0.2);
    p3->SetLogz();
    p3->Draw();
    p3->cd();
    h_response->SetTitle("");
    h_response->SetStats(0);
    h_response->GetXaxis()->SetTitle("Reconstructed jet #font[52]{p}_{T} [GeV]");
    h_response->GetYaxis()->SetTitle("Truth jet #font[52]{p}_{T} [GeV]");
    h_response->Draw("colz");
    h_response->GetZaxis()->SetRangeUser(0.000001,7e5);
   


    









}
