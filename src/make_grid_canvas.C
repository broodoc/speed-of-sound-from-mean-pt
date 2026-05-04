#include <TCanvas.h>
#include <TPad.h>
#include <TImage.h>
#include <TString.h>
#include <vector>

void make_grid_canvas(const char* imageDir = "./") {
    const int nCols = 6;
    const int nRows = 6;

    std::vector<TString> images = {
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.000-0.001%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.001-0.002%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.002-0.003%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.003-0.005%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.005-0.010%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.010-0.020%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.020-0.030%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.030-0.050%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.050-0.100%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.100-0.200%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.200-0.300%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.300-0.500%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__0.500-1.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__1.000-2.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__2.000-3.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__3.000-4.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__4.000-5.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__5.000-5.500%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__5.500-6.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__6.000-7.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__7.000-8.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__8.000-9.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__9.000-10.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__10.000-15.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__15.000-20.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__20.000-25.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__25.000-30.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__30.000-40.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__40.000-50.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__50.000-60.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__60.000-70.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__70.000-80.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__80.000-90.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__90.000-100.000%_LevyTsallis.png",
        "BestChi2Ndf_spectra_study_fit_range_0.00-1.25_GeV__ref_0-5%_LevyTsallis.png"
    };

    const int nPerPage = nCols * nRows;
    const int nPages = (images.size() + nPerPage - 1) / nPerPage;

    for (int page = 0; page < nPages; ++page) {
        TCanvas* c = new TCanvas(Form("c_page_%d", page), "Grid", 2400, 1600);
        c->Divide(nCols, nRows);

        for (int i = 0; i < nPerPage; ++i) {
            int idx = page * nPerPage + i;
            if (idx >= images.size()) break;

            TString path = Form("%s/%s", imageDir, images[idx].Data());
            TImage* img = TImage::Open(path);
            if (!img) {
                printf("Failed to open %s\n", path.Data());
                continue;
            }

            c->cd(i + 1);
            img->Draw("xxx");  // draw into pad
            delete img;
        }

        c->SaveAs(Form("grid_canvas_page_%02d.pdf", page + 1));
        delete c;
    }
}
