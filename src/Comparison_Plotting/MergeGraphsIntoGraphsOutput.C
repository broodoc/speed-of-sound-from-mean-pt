#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>

void MergeGraphsIntoGraphsOutput() {
    // Open ROOT files
    TFile *f1 = TFile::Open("/Users/calebbroodo/Desktop/sos_analysis/src/ana_code/graphs_output.root", "UPDATE"); // Open in update mode
//    TFile *f2 = TFile::Open("Trajectum_output.root", "READ");
    TFile *f3 = TFile::Open("CMS_SOS.root", "READ");
    TFile *f4 = TFile::Open("HIJING_Output.root", "READ");
    TFile *f5 = TFile::Open("/Users/calebbroodo/Desktop/sos_analysis/prepped_files/54p5COL_Run18/54p5COL_Run18_refMult2/graphs_output_2025-Apr-29_11-14AM.root", "READ");
//    TFile *f6 = TFile::Open("/Users/calebbroodo/sos_analysis/prepped_files/200COL_Run10_Run11/200COL_Run11_Run10_refMult2/graphs_output_2025-Apr-28_12-41PM.root", "READ");
    TFile *f7 = TFile::Open("Trajectum_output_absY_normX.root", "READ");
    TFile *f8 = TFile::Open("PHOBOS_200GeV_Measurement_arXiv_nucl_ex_0201005v2.root", "READ");
    TFile *f9 = TFile::Open("STAR_200GeV_Measurement_doi_10.1088slash1742-6596slash5slash1slash003.root", "READ");
    TFile *f10 = TFile::Open("Traj_200GeV_Measurement_July17_centrality_pTlow0p20_InclusivepT.root", "READ");
    

//    if (!f1 || !f3 || !f4 ||!f5 || !f6 || !f7 || !f8 || !f9 || !f10) {
//        std::cerr << "❌ Error: Could not open one or more ROOT files!" << std::endl;
//        return;
//    }
    if (!f1 || !f3 || !f4 ||!f5 || !f7 || !f8 || !f9 || !f10) {
        std::cerr << "❌ Error: Could not open one or more ROOT files!" << std::endl;
        return;
    }

    std::cout << "📈 Retrieving graphs for comparison" << std::endl;
    // Retrieve graphs from Trajectum and CMS files
//    TGraph *g2 = (TGraph*)f2->Get("Mean_pT_vs_Centrality");
    TGraph *g3 = (TGraph*)f3->Get("Figure 3/Graph1D_y1");
    TGraph *g4 = (TGraph*)f4->Get("HIJING Au+Au #sqrt{s_{NN}} = 200 GeV (y2016)");
    TGraph *g7 = (TGraph*)f7->Get("Mean_pT_vs_Centrality");
    TGraph *g8 = (TGraph*)f8->Get("h_dNdeta");
    TGraph *g9 = (TGraph*)f9->Get("h_dNdeta");
    TGraph *g10 = (TGraph*)f10->Get("g_dNdeta");
    std::cout << "📈 Retrieving 54.4 GeV for comparison" << std::endl;
    TDirectory *gdir = (TDirectory*)f5->Get("pT study/Normalized (to ref class) quantities/pT-Integrated Nch/fit_range_0.00-1.45_GeV");
    TKey *key = gdir->GetKey("0 #leq p_{T} #leq 10 (GeV/c)");
    TGraph *g5 = (TGraphErrors*)key->ReadObject<TGraphErrors>();
    
    std::cout << "📈 Retrieving 200 GeV refMult2 for comparison" << std::endl;
//    TDirectory *gdir1 = (TDirectory*)f6->Get("pT study/Normalized (to ref class) quantities/pT-Integrated Nch/fit_range_0.00-0.40_GeV");
//    TKey *key1 = gdir1->GetKey("0 #leq p_{T} #leq 10 (GeV/c)");
//    TGraph *g6 = (TGraphErrors*)key1->ReadObject<TGraphErrors>();
    
//    TDirectory *gdir2 = (TDirectory*)f6->Get("General");
//    TKey *key2 = gdir2->GetKey("Centrality__pTlow0p20_hRefMult");
//    TGraph *g7 = (TGraphErrors*)key2->ReadObject<TGraphErrors>();
    
    std::cout << "📈 Cloning and writing graphs to graphs_output.root" << std::endl;
//    if (!g2) std::cerr << "❌ Error: Could not retrieve Trajectum graph!" << std::endl;
    if (!g3) std::cerr << "❌ Error: Could not retrieve CMS graph!" << std::endl;
    if (!g4) std::cerr << "❌ Error: Could not retrieve HIJING graph!" << std::endl;
    if (!g5) std::cerr << "❌ Error: Could not retrieve 54.4 GeV refMult2 graph!" << std::endl;
//    if (!g6) std::cerr << "❌ Error: Could not retrieve 200 GeV refMult2 graph!" << std::endl;
    if (!g7) std::cerr << "❌ Error: Could not retrieve Trajectum (absY) graph!" << std::endl;
//    if (!g7) std::cerr << "❌ Error: Could not retrieve 200 GeV refMult2 Nch graph!" << std::endl;

    // Create a directory for new graphs if it doesn't exist
    TDirectory *mergedDir = f1->GetDirectory("Merged_Graphs");
    if (!mergedDir) {
        mergedDir = f1->mkdir("Merged_Graphs");
    }
    mergedDir->cd(); // Move into the directory

    // Clone and write graphs if they exist
//    if (g2) {
//        TGraph *g2_clone = (TGraph*)g2->Clone("Trajectum_Graph");
//        g2_clone->SetTitle("Trajectum (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} #leq 10 (GeV/c), |#eta| < 0.5");
//        g2_clone->SetLineWidth(2);
//        g2_clone->SetLineColor(kBlack);
//        g2_clone->SetMarkerStyle(20);
//        g2_clone->SetMarkerColor(kAzure+1);
//        g2_clone->Write();
//        std::cout << "📈 Added\033[34m Trajectum\033[0m graph to graphs_output.root" << std::endl;
//    }

    if (g3) {
        TGraph *g3_clone = (TGraph*)g3->Clone("CMS_Graph");
        g3_clone->SetTitle("CMS (Pb+Pb #sqrt{s_{NN}} = 5.02 TeV), p_{T} > 0 (GeV/c), |#eta| < 0.5");
        g3_clone->SetLineWidth(2);
        g3_clone->SetLineColor(kBlack);
        g3_clone->SetMarkerStyle(20);
        g3_clone->SetMarkerColor(kRed);
        g3_clone->Write();
        std::cout << "📈 Added\033[31m CMS\033[0m graph to graphs_output.root" << std::endl;
    }
    if (g4) {
        TGraph *g4_clone = (TGraph*)g4->Clone("HIJING_Graph");
        g4_clone->SetTitle("HIJING (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} > 0 (GeV/c), |#eta| < 0.5");
        g4_clone->SetLineWidth(2);
        g4_clone->SetLineColor(kBlack);
        g4_clone->SetMarkerStyle(20);
        g4_clone->SetMarkerColor(kMagenta);
        g4_clone->Write();
        std::cout << "📈 Added\033[35m HIJING\033[0m graph to graphs_output.root" << std::endl;
    }
    if (g5) {
        TGraph *g5_clone = (TGraph*)g5->Clone("refMult2_Graph 52p4 GeV");
        g5_clone->SetTitle("refMult2 (Au+Au #sqrt{s_{NN}} = 54.4 GeV), p_{T} > 0 (GeV/c), |#eta| < 0.5");
        g5_clone->SetLineWidth(2);
        g5_clone->SetLineColor(kBlack);
        g5_clone->SetMarkerStyle(20);
        g5_clone->SetMarkerColor(kMagenta);
        g5_clone->Write();
        std::cout << "📈 Added\033[35m refMult2\033[0m graph to graphs_output.root" << std::endl;
    }
//    if (g6) {
//        TGraph *g6_clone = (TGraph*)g6->Clone("refMult2_Graph 200 GeV");
//        g6_clone->SetTitle("refMult2 (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} > 0 (GeV/c), |#eta| < 0.5");
//        g6_clone->SetLineWidth(2);
//        g6_clone->SetLineColor(kBlack);
//        g6_clone->SetMarkerStyle(20);
//        g6_clone->SetMarkerColor(kMagenta);
//        g6_clone->Write();
//        std::cout << "📈 Added\033[35m refMult2\033[0m graph to graphs_output.root" << std::endl;
//    }
    if (g7) {
        TGraph *g7_clone = (TGraph*)g7->Clone("Trajectum_absY_normX_Graph");
        g7_clone->SetTitle("Trajectum (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} #leq 10 (GeV/c), |#eta| < 0.5");
        g7_clone->SetLineWidth(2);
        g7_clone->SetLineColor(kBlack);
        g7_clone->SetMarkerStyle(20);
        g7_clone->SetMarkerColor(kAzure+1);
        g7_clone->Write();
        std::cout << "📈 Added\033[34m Trajectum (absY)\033[0m graph to graphs_output.root" << std::endl;
    }
    if (g8) {
        TGraph *g8_clone = (TGraph*)g8->Clone("PHOBOS 2002 Au+Au 200 GeV dN/d#eta");
        g8_clone->SetTitle("PHOBOS (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} #leq 10 (GeV/c), |#eta| < 0.5");
        g8_clone->SetLineColor(kRed);
        g8_clone->SetMarkerStyle(20);
        g8_clone->SetMarkerColor(kRed+1);
        g8_clone->Write("PHOBOS 2002 Au+Au 200 GeV dN/d#eta");
        std::cout << "📈 Added\033[34m PHOBOS\033[0m graph to graphs_output.root" << std::endl;
        }
    if (g9) {
        TGraph *g9_clone = (TGraph*)g9->Clone("STAR Au+Au 200 GeV 2005 <pT> dN/d#eta midrapidity");
        g9_clone->SetTitle("STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} #leq 10 (GeV/c), |#eta| < 0.5");
        g9_clone->SetLineColor(kRed);
        g9_clone->SetMarkerStyle(20);
        g9_clone->SetMarkerColor(kRed+1);
        g9_clone->Write("STAR Au+Au 200 GeV 2005 <pT> dN/d#eta midrapidity");
        std::cout << "📈 Added\033[34m STAR\033[0m graph to graphs_output.root" << std::endl;
        }
    if (g10 && g10->GetN() > 0) {
        // Clone the graph
        TGraph* g10_clone = (TGraph*)g10->Clone("Trajectum Au+Au 200 GeV <pT> dN/d#eta midrapidity (July 17)");

        // Define the PowerLaw fit function
        TF1* powerLawFit = new TF1("PowerLaw_fit", "[1]*pow(x, [0])", 1.14, 1.3);
        powerLawFit->SetParameter(0, 0.2);  // Initial guess for exponent (cs^2)
        powerLawFit->SetParameter(1, 1.0);  // Initial guess for normalization

        // Fit with fit display (option "SMREQ" = Show, Minos errors, Quiet)
        g10_clone->Fit(powerLawFit, "SMREQ");

        // Extract fit results
        double cs2     = powerLawFit->GetParameter(0);
        double cs2Err  = powerLawFit->GetParError(0);
        double chi2ndf = (powerLawFit->GetNDF() > 0) ? powerLawFit->GetChisquare() / powerLawFit->GetNDF() : -1;

        // Create a legend
        TLegend* leg = new TLegend(0.15, 0.65, 0.6, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->AddEntry(g10_clone, "STAR Au+Au  #sqrt{s_{NN}} = 200 GeV", "p");
        leg->AddEntry((TObject*)0, "0.0 #leq p_{T} #leq 10.0 GeV/c, |#eta| #leq 0.5", "");
        leg->AddEntry((TObject*)0, Form("#font[12]{Fit to extract}  c_{s}^{2} [1.075, 1.250]"), "");
        leg->AddEntry((TObject*)0, Form("c_{s}^{2} = %.4f #pm %.6f (stat)", cs2, cs2Err), "");
        leg->AddEntry((TObject*)0, Form("#chi^{2} / NDF = %.2f", chi2ndf), "");

        // Set style
        g10_clone->SetTitle("Trajectum (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} #leq 10 (GeV/c), |#eta| < 0.5");
        g10_clone->SetLineColor(kBlue+1);
        g10_clone->SetMarkerStyle(20);
        g10_clone->SetMarkerColor(kBlue+1);

        // Draw everything to a canvas
        TCanvas* c1 = new TCanvas("c1", "Fit with annotation", 800, 600);
        g10_clone->Draw("AP");
        powerLawFit->SetLineColor(kBlack);
        powerLawFit->Draw("same");
        leg->Draw();

        // Save to output
        g10_clone->Write("Trajectum Au+Au 200 GeV <pT> dN/d#eta midrapidity (July 17)");
        c1->Write("FitWithLegend");

        std::cout << "🔧 Fit: c_s^2 = " << cs2 << " ± " << cs2Err << ", χ²/NDF = " << chi2ndf << std::endl;
        std::cout << "📈 Added\033[34m Trajectum\033[0m graph + fit to graphs_output.root" << std::endl;
    }

    
//    if (g10) {
//        TGraph *g10_clone = (TGraph*)g10->Clone("Trajectum Au+Au 200 GeV <pT> dN/d#eta midrapidity (July 17)");
//        g10_clone->SetTitle("Trajectum (Au+Au #sqrt{s_{NN}} = 200 GeV), p_{T} #leq 10 (GeV/c), |#eta| < 0.5");
//        g10_clone->SetLineColor(kRed);
//        g10_clone->SetMarkerStyle(20);
//        g10_clone->SetMarkerColor(kBlue+1);
//        g10_clone->Write("Trajectum Au+Au 200 GeV <pT> dN/d#eta midrapidity (July 17)");
//        std::cout << "📈 Added\033[34m Trajectum\033[0m graph to graphs_output.root" << std::endl;
//        }
    
    // Finalize
    f1->Write();
    f1->Close();
//    f2->Close();
    f3->Close();
    f4->Close();
    f5->Close();
//    f6->Close();
    f8->Close();
    f9->Close();
    f10->Close();
    std::cout << "✅ Successfully merged graphs into graphs_output.root" << std::endl;
}



