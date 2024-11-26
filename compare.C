/*
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

void compare() {
  TFile *file1 = TFile::Open("/home/shalini/finaldestiny/atlaspaper_0.root");
  TFile *file2 = TFile::Open("/home/shalini/finaldestiny/atlaspaper_1.root");
  TFile *file3 = TFile::Open("/home/shalini/finaldestiny/atlaspaper_27.root");
  
  // Retrieve aplanarityHist histograms from each file
  TH1F *aplanarityHist1 = (TH1F*)file1->Get("aplanarityHist");
  TH1F *aplanarityHist2 = (TH1F*)file2->Get("aplanarityHist");
  TH1F *aplanarityHist3 = (TH1F*)file3->Get("aplanarityHist");

  TH1F *sphericityHist1 = (TH1F*)file1->Get("sphericityHist");
  TH1F *sphericityHist2 = (TH1F*)file2->Get("sphericityHist");
  TH1F *sphericityHist3 = (TH1F*)file3->Get("sphericityHist");

  
  // Clone histograms to avoid modifying originals
  TH1F *histAplanarity1 = (TH1F*)aplanarityHist1->Clone("histAplanarity1");
  TH1F *histAplanarity2 = (TH1F*)aplanarityHist2->Clone("histAplanarity2");
  TH1F *histAplanarity3 = (TH1F*)aplanarityHist3->Clone("histAplanarity3");


  TH1F *histSphericity1 = (TH1F*)sphericityHist1->Clone("histSphericity1");
  TH1F *histSphericity2 = (TH1F*)sphericityHist2->Clone("histSphericity2");
  TH1F *histSphericity3 = (TH1F*)sphericityHist3->Clone("histSphericity3");
  
 // Book the histograms (aplanarity)
    histAplanarity1->SetTitle("Aplanarity Comparison; Aplanarity; Events");
    histAplanarity1->SetLineColor(kRed);
    histAplanarity1->SetMarkerColor(kRed);
    histAplanarity1->SetMarkerStyle(20);  // Circle marker
    histAplanarity1->SetMarkerSize(1);

    histAplanarity2->SetLineColor(kBlue);
    histAplanarity2->SetMarkerColor(kBlue);
    histAplanarity2->SetMarkerStyle(21);  // Square marker
    histAplanarity2->SetMarkerSize(1);

    histAplanarity3->SetLineColor(kGreen);
    histAplanarity3->SetMarkerColor(kGreen);
    histAplanarity3->SetMarkerStyle(22);  // Triangle marker
    histAplanarity3->SetMarkerSize(1);

    // Book the histograms (sphericity)
    histSphericity1->SetTitle("Sphericity Comparison; Sphericity; Events");
    histSphericity1->SetLineColor(kRed);
    histSphericity1->SetMarkerColor(kRed);
    histSphericity1->SetMarkerStyle(24);  // Cross marker
    histSphericity1->SetMarkerSize(1);

    histSphericity2->SetLineColor(kBlue);
    histSphericity2->SetMarkerColor(kBlue);
    histSphericity2->SetMarkerStyle(25);  // Circle marker
    histSphericity2->SetMarkerSize(1);

    histSphericity3->SetLineColor(kGreen);
    histSphericity3->SetMarkerColor(kGreen);
    histSphericity3->SetMarkerStyle(26);  // Triangle marker
    histSphericity3->SetMarkerSize(1);

    // Create canvas with 2 pads for separate plots
    TCanvas *c1 = new TCanvas("c1", "Comparison of Aplanarity and Sphericity", 1200, 800);
    c1->Divide(1, 2);  // Divide canvas into 2 parts (1 row, 2 columns)

    // First pad: Aplanarity
    c1->cd(1);  // Switch to the first pad
    histAplanarity1->Draw("E");
    histAplanarity2->Draw("E SAME");
    histAplanarity3->Draw("E SAME");

    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend1->AddEntry(histAplanarity1, "Energy 1", "p");
    legend1->AddEntry(histAplanarity2, "Energy 2", "p");
    legend1->AddEntry(histAplanarity3, "Energy 3", "p");
    legend1->Draw();

    // Second pad: Sphericity
    c1->cd(2);  // Switch to the second pad
    histSphericity1->Draw("E");
    histSphericity2->Draw("E SAME");
    histSphericity3->Draw("E SAME");

    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(histSphericity1, "Energy 1", "p");
    legend2->AddEntry(histSphericity2, "Energy 2", "p");
    legend2->AddEntry(histSphericity3, "Energy 3", "p");
    legend2->Draw();

    // Update canvas
    c1->Update();
}


*/

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

void compare() {
  // Open the ROOT files
  TFile *file1 = TFile::Open("/home/shalini/finaldestiny/atlaspaper_27mpioff.root");
  TFile *file2 = TFile::Open("/home/shalini/finaldestiny/atlaspaper_13.6mpioff.root");
  TFile *file3 = TFile::Open("/home/shalini/finaldestiny/atlaspaper_8mpioff.root");

  // Retrieve histograms from each file
  TH1F *aplanarityHist1 = (TH1F*)file1->Get("aplanarityHist");
  TH1F *aplanarityHist2 = (TH1F*)file2->Get("aplanarityHist");
  TH1F *aplanarityHist3 = (TH1F*)file3->Get("aplanarityHist");

  TH1F *sphericityHist1 = (TH1F*)file1->Get("sphericityHist");
  TH1F *sphericityHist2 = (TH1F*)file2->Get("sphericityHist");
  TH1F *sphericityHist3 = (TH1F*)file3->Get("sphericityHist");

  TH1F *transSphericityHist1 = (TH1F*)file1->Get("transSphericityrHist");
  TH1F *transSphericityHist2 = (TH1F*)file2->Get("transSphericityrHist");
  TH1F *transSphericityHist3 = (TH1F*)file3->Get("transSphericityrHist");

  // Normalize the histograms to compare shapes
  aplanarityHist1->Scale(1.0 / aplanarityHist1->Integral());
  aplanarityHist2->Scale(1.0 / aplanarityHist2->Integral());
  aplanarityHist3->Scale(1.0 / aplanarityHist3->Integral());

  sphericityHist1->Scale(1.0 / sphericityHist1->Integral());
  sphericityHist2->Scale(1.0 / sphericityHist2->Integral());
  sphericityHist3->Scale(1.0 / sphericityHist3->Integral());

  transSphericityHist1->Scale(1.0 / transSphericityHist1->Integral());
  transSphericityHist2->Scale(1.0 / transSphericityHist2->Integral());
  transSphericityHist3->Scale(1.0 / transSphericityHist3->Integral());

  // Set line colors and styles
  aplanarityHist1->SetLineColor(kRed);
  aplanarityHist2->SetLineColor(kBlue);
  aplanarityHist3->SetLineColor(kGreen);

  sphericityHist1->SetLineColor(kRed);
  sphericityHist2->SetLineColor(kBlue);
  sphericityHist3->SetLineColor(kGreen);

  transSphericityHist1->SetLineColor(kRed);
  transSphericityHist2->SetLineColor(kBlue);
  transSphericityHist3->SetLineColor(kGreen);

  // Create canvases for each histogram
  TCanvas *c1 = new TCanvas("c1", "Comparison of Aplanarity", 1200, 800);
  aplanarityHist1->Draw("HIST");
  aplanarityHist2->Draw("HIST SAME");
  aplanarityHist3->Draw("HIST SAME");
  
  // Add legend for aplanarity plot
  TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->AddEntry(aplanarityHist1, "27 TeV", "l");
  legend1->AddEntry(aplanarityHist2, "13.6 TeV", "l");
  legend1->AddEntry(aplanarityHist3, "8 TeV", "l");
  legend1->Draw();
  c1->SetLogy();  // Set logarithmic y-axis
  c1->Update();

  // Create canvas for Sphericity
  TCanvas *c2 = new TCanvas("c2", "Comparison of Sphericity", 1200, 800);
  sphericityHist1->Draw("HIST");
  sphericityHist2->Draw("HIST SAME");
  sphericityHist3->Draw("HIST SAME");
  
  // Add legend for sphericity plot
  TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2->AddEntry(sphericityHist1, "27 TeV", "l");
  legend2->AddEntry(sphericityHist2, "13.6 TeV", "l");
  legend2->AddEntry(sphericityHist3, "8 TeV", "l");
  legend2->Draw();
  c2->SetLogy();  // Set logarithmic y-axis
  c2->Update();

  // Create canvas for Transverse Sphericity
  TCanvas *c3 = new TCanvas("c3", "Comparison of Transverse Sphericity", 1200, 800);
  transSphericityHist1->Draw("HIST");
  transSphericityHist2->Draw("HIST SAME");
  transSphericityHist3->Draw("HIST SAME");
  
  // Add legend for transverse sphericity plot
  TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend3->AddEntry(transSphericityHist1, "27 TeV", "l");
  legend3->AddEntry(transSphericityHist2, "13.6 TeV", "l");
  legend3->AddEntry(transSphericityHist3, "8 TeV", "l");
  legend3->Draw();
  c3->SetLogy();  // Set logarithmic y-axis
  c3->Update();
}
