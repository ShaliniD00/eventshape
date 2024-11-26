/*
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

using namespace std;

void version8oct() {
    // Open the ROOT file
    TFile *file = new TFile("/home/shalini/finaldestiny/jetplots_13mpion_improved.root", "READ");
    if (!file->IsOpen()) {
        cout << "Error opening file" << endl;
        return;
    }

    // List of histogram names
    const char* histNames[] = {
        "sphericityHist_AllJets",
        "sphericityHist_0Jets",
        "sphericityHist_1Jet",
        "sphericityHist_2Jets",
        "sphericityHist_3PlusJets"
    };

    // Colors for each histogram
    Color_t colors[] = {kBlack, kRed, kBlue, kGreen, kMagenta};

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Overlaid Histograms", 800, 600);

    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Create a histogram stack
    THStack *stack = new THStack("stack", "Transverse Sphericity (All Jets)");

    double maxY = 0;

    // Get and prepare histograms
    for (int i = 0; i < 5; ++i) {
        TH1F *hist = (TH1F*)file->Get(histNames[i]);
        if (!hist) {
            cout << "Histogram " << histNames[i] << " not found in the file." << endl;
            continue;
        }

        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(2);

        // Find the maximum y-value across all histograms
        double histMax = hist->GetMaximum();
        if (histMax > maxY) maxY = histMax;

        stack->Add(hist);

        TString legendEntry = histNames[i];
        legendEntry = legendEntry(legendEntry.Last('_')+1, legendEntry.Length());
        legend->AddEntry(hist, legendEntry, "l");
    }

    // Draw the stack
    stack->Draw("nostack");
    stack->SetMaximum(maxY * 1.1);  // Set the maximum y-value with some headroom

    // Set axis labels
    stack->GetXaxis()->SetTitle("S_{T}");
    stack->GetYaxis()->SetTitle("dN/dS_{T}");

    // Set y-axis to log scale
    canvas->SetLogy();

    // Draw the legend
    legend->Draw();

    // Update the canvas
    canvas->Update();

    // Save the plot
    canvas->SaveAs("overlaid_sphericity_histograms.png");

    // Keep the canvas open (for interactive sessions)
    canvas->Draw();

    // Close the file
    //file->Close();
}

*/

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <THStack.h>
#include <iostream>

using namespace std;

void version8oct() {
    // Open the ROOT file
    TFile *file = new TFile("/home/shalini/finaldestiny/jetplots_13mpion_improved.root", "READ");
    if (!file->IsOpen()) {
        cout << "Error opening file" << endl;
        return;
    }

    // List of histogram names
    const char* histNames[] = {
        "sphericityHist_AllJets",
        "sphericityHist_0Jets",
        "sphericityHist_1Jet",
        "sphericityHist_2Jets",
        "sphericityHist_3PlusJets"
    };

    // Colors for each histogram
    Color_t colors[] = {kBlack, kRed, kBlue, kGreen, kMagenta};

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Overlaid Histograms", 800, 600);

    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Create a histogram stack
    THStack *stack = new THStack("stack", "Transverse Sphericity (All Jets)");

    // Get and prepare histograms
    for (int i = 0; i < 5; ++i) {
        TH1F *hist = (TH1F*)file->Get(histNames[i]);
        if (!hist) {
            cout << "Histogram " << histNames[i] << " not found in the file." << endl;
            continue;
        }

        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(2);
        
        // Find the maximum y-value across all histograms
        stack->Add(hist);

        TString legendEntry = histNames[i];
        legendEntry = legendEntry(legendEntry.Last('_') + 1, legendEntry.Length());
        legend->AddEntry(hist, legendEntry, "l");
    }

    // Draw the stack
    stack->Draw("nostack");

    // Set axis labels
    stack->GetXaxis()->SetTitle("S_{T}");
    stack->GetYaxis()->SetTitle("dN/dS_{T}");

    // Set y-axis to log scale
    canvas->SetLogy();

    // Draw the 0-jet histogram separately to set its range
    TH1F *zeroJetHist = (TH1F*)file->Get(histNames[1]);
    if (zeroJetHist) {
        zeroJetHist->SetLineColor(colors[1]); // Ensure the color matches
        zeroJetHist->SetLineWidth(2);
        zeroJetHist->Draw("same"); // Draw it on top of the stack
    }

    // Set the y-axis limits specifically for the 0-jet histogram
    double yMin = 1e-3;
    double yMax = 1e-1;

    // Get the current maximum y value of the stack
    double maxY = stack->GetMaximum();
    if (maxY > yMax) {
        stack->SetMaximum(maxY * 1.1); // Set maximum slightly above the max value of stack
    } else {
        stack->SetMaximum(yMax); // Otherwise, set to the specified maximum
    }
    stack->GetYaxis()->SetRangeUser(yMin, yMax); // Set the y-axis range

    // Draw the legend
    legend->Draw();

    // Update the canvas
    canvas->Update();

    // Save the plot
    canvas->SaveAs("overlaid_sphericity_histograms.png");

    // Keep the canvas open (for interactive sessions)
    canvas->Draw();

    // Close the file
    // file->Close();
}
