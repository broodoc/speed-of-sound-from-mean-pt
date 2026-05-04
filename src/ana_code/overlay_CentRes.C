#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>

void overlay() {
    gROOT->SetBatch(kTRUE);
    TFile* file = TFile::Open("graphs_dEta1p0__pTlow0p20_hPtDec_Run21_from_raw_spectrum.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return;
    }

    TH1D* h_black = (TH1D*)file->Get("h_black");
    TH1D* h_red = (TH1D*)file->Get("Analytic Model: c_{s}^{2} = 0.08");
    TGraphErrors* star = (TGraphErrors*)file->Get("STAR Au+Au #sqrt{s_{NN}} = 200 GeV Run16");
    TGraphErrors* c_r = (TGraphErrors*)file->Get("centrality_resolution_Nch");
    

    if (!h_black || !c_r) {
        std::cerr << "Error: Could not retrieve h_black or c_r." << std::endl;
        file->Close();
        return;
    }
    
    for (int i = 1; i <= h_black->GetNbinsX(); ++i) {
        if (i % 1000 == 0) std::cout << "Bin " << i << ": Content = " << h_black->GetBinContent(i) << std::endl;
    }


    TCanvas* c = new TCanvas("c", "Overlay with Dual Y-Axis", 800, 600);
    c->SetCanvasSize(800, 600);


    h_black->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_black->GetYaxis()->SetTitle("");
    h_black->GetXaxis()->SetTitle("N_{ch}/N_{ch,0-5%}");
    h_black->SetTitle("");
    h_black->SetLineColor(kBlack);
    h_black->SetLineWidth(5);
    h_black->SetStats(0);
    h_black->Draw("HIST");
    gPad->Modified();
    gPad->Update();
    
    h_red->GetYaxis()->SetRangeUser(0.8, 1.1);
    h_red->GetYaxis()->SetTitle("");
    h_red->GetXaxis()->SetTitle("N_{ch}/N_{ch,0-5%}");
    h_red->SetTitle("");
    h_red->SetLineColor(kRed);
    h_red->SetLineWidth(5);
    h_red->SetStats(0);
    h_red->Draw("HIST SAME");
    

//    star->SetLineColor(kRed);
//    star->SetLineWidth(5);
    star->SetStats(0);
    star->Draw("HIST SAME");
    
    TLegend* legend = new TLegend(0.2, 0.7, 0.5, 0.8); // Define legend position (top left corner)
    legend->SetBorderSize(0); // Remove the outline
    legend->SetFillStyle(0);  // Make the background transparent
    legend->SetTextSize(0.03); // Reduce the font size
    legend->AddEntry(h_black, "f(N_{ch}/N_{ch,0-5%})", "l"); // Label for h_black
    legend->AddEntry(h_red, "#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}", "l"); // Label for h_red
    legend->AddEntry(star, "STAR Au+Au #sqrt{s_{NN}} = 200 GeV Run16", "l"); // Label for star
    legend->Draw();
    
    gPad->Modified();
    gPad->Update();

    c->Update();
    c->Draw();

    //    double scale = (1.1 - 0.8) / (10.0 - 0.0);
        
    int nPoints = c_r->GetN();
    double y_min_left = 0.8;
    double y_range_left = 0.3;
    double y_min_right = 0;
    double y_range_right = 10;
    double scale_factor = y_range_left / y_range_right;
    double y_line = 5 * scale_factor + y_min_left;
    
        TLine* hLine = new TLine(h_black->GetXaxis()->GetXmin(), y_line,
                                 h_black->GetXaxis()->GetXmax(), y_line);
        hLine->SetLineStyle(2); // Dashed line
        hLine->SetLineColor(kOrange - 3); // Choose a color for visibility
        hLine->Draw();
    
        TLine* vLine = new TLine(1.1156, gPad->GetUymin(),
                                 1.1156, gPad->GetUymax());
        vLine->SetLineStyle(2); // Dashed line
        vLine->SetLineColor(kBlack); // Choose a color for visibility
        vLine->Draw();
    
    TLatex* label = new TLatex(1.113, (gPad->GetUymax() - gPad->GetUymin()) / 1.75 + gPad->GetUymin(), "knee");
    label->SetTextSize(0.03);
    label->SetTextColor(kBlack);
    label->SetTextAngle(90);
    label->SetTextFont(42);
    label->Draw();
    
    for (int i = 0; i < nPoints; ++i) {
        double x, y, ex, ey;
        c_r->GetPoint(i, x, y);

        // Scale the y value
        double y_temp = y; // Store original y value
        y = y_min_left + (y - y_min_right) * (y_range_left / y_range_right);

        c_r->SetPoint(i, x, y);
//        ex = 0.1;
//        ey = 0;
//
//        c_r->SetPointError(i, ex, ey);
        cout << "i: " << i << " | x: " << x << " | y (scaled): " << y
             << " | y (original): " << y_temp << endl;
    }

    
    double right_y_min = 0;
    double right_y_max = 10;

    TGaxis* rightAxis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                   gPad->GetUxmax(), gPad->GetUymax(),
                                   right_y_min, right_y_max, 510, "+L");
//    rightAxis->SetLineColor(kGray + 2);
    rightAxis->SetLabelColor(kOrange - 3);
    rightAxis->SetTitleColor(kOrange - 3);
    rightAxis->SetTitle("Relative centrality error (%)");
    rightAxis->Draw();
    
    TGaxis* upperAxis = new TGaxis(0.8, gPad->GetUymax(), 1.3, gPad->GetUymax(), 0.8, 1.3, 510, "-U");
        upperAxis->SetTitleOffset(1.2);
        upperAxis->Draw();

//    c_r->SetLineColor(kMagenta);      // Set the line color to black
//    c_r->SetFillColorAlpha(kMagenta, 0.22); // Set the fill color to black with opacity 0.22
//    c_r->SetMarkerStyle(20);        // Use a marker style for additional visibility (optional)
//    c_r->Draw("B SAME");
    
    TGraphErrors* c_r_points = (TGraphErrors*)c_r->Clone(); // Clone to preserve properties
    c_r_points->SetMarkerStyle(20);  // Set marker style
    c_r_points->SetMarkerSize(1.5);
    c_r_points->SetMarkerColor(kOrange - 3); // Set marker color (e.g., red for visibility)
    c_r_points->Draw("P SAME");       // Draw only the points
    
//    c_r->SetMarkerStyle(20);
//    c_r->SetMarkerColor(kGray + 2);
//    c_r->SetLineColor(kGray + 2);
//    c_r->Draw("P SAME");
    

//    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
//    legend->AddEntry(h_black, "f(N_{ch}/N_{ch,0-5%})", "l");
//    legend->AddEntry(c_r, "Centrality Resolution", "p");
//    legend->Draw();

    std::cout << "Left axis range: " << h_black->GetMinimum() << ", " << h_black->GetMaximum() << std::endl;
    std::cout << "Right axis range: " << right_y_min << ", " << right_y_max << std::endl;

    c->SaveAs("Graphs/overlay_with_star_black_red_and_orange_fixed.pdf");
    file->Close();
}
