#ifndef MOMENTHISTOGRAM_H
#define MOMENTHISTOGRAM_H

#include "Parameterization.h"

class MomentHistogram : public TObject {
    
public:
    MomentHistogram();
    ~MomentHistogram();
    
    void InitHistograms();
    
    //Curves and calculations guided by page 9 of
    //https://arxiv.org/pdf/2101.12413
    
    // -----**-----**-----**-----**-----**--p_{T} quantities--**-----**-----**-----**-----**
    
    double pTC1_ref;
    double pTC2_ref;
    double pTC3_ref;
    double pTC4_ref;
    
    double pTC1;
    double pTC2;
    double pTC3;
    double pTC4;

    double pTM_ref;
    double pTVar_ref;
    double pTS_ref;
    double pTkap_ref;
    
    double pTM;
    double pTVar;
    double pTS;
    double pTkap;
    
    // -----**-----**-----**-----**-----**--N_{ch} quantities--**-----**-----**-----**-----**
    
    double NchC1_ref;
    double NchC2_ref;
    double NchC3_ref;
    double NchC4_ref;

    double NchC1;
    double NchC2;
    double NchC3;
    double NchC4;

    double NchM_ref;
    double NchVar_ref;
    double NchS_ref;
    double Nchkap_ref;

    double NchM;
    double NchVar;
    double NchS;
    double Nchkap;
    
    // -----**-----**-----**-----**-----**--p_{T} quantities--**-----**-----**-----**-----**

    double pTC1_ref_error;
    double pTC2_ref_error;
    double pTC3_ref_error;
    double pTC4_ref_error;

    double pTC1_error;
    double pTC2_error;
    double pTC3_error;
    double pTC4_error;

    double pTM_ref_error;
    double pTVar_ref_error;
    double pTS_ref_error;
    double pTkap_ref_error;

    double pTM_error;
    double pTVar_error;
    double pTS_error;
    double pTkap_error;

    // -----**-----**-----**-----**-----**--N_{ch} quantities--**-----**-----**-----**-----**

    double NchC1_ref_error;
    double NchC2_ref_error;
    double NchC3_ref_error;
    double NchC4_ref_error;

    double NchC1_error;
    double NchC2_error;
    double NchC3_error;
    double NchC4_error;

    double NchM_ref_error;
    double NchVar_ref_error;
    double NchS_ref_error;
    double Nchkap_ref_error;

    double NchM_error;
    double NchVar_error;
    double NchS_error;
    double Nchkap_error;


    // -----**-----**-----**-----**-----**--p_{T} Graphs--**-----**-----**-----**-----**
    
    // Cumulants (Absolute)
    TGraphErrors* gpTC1;
    TGraphErrors* gpTC2;
    TGraphErrors* gpTC3;
    TGraphErrors* gpTC4;
    
    // Cumulants (Normalized to 0-5%)
    TGraphErrors* gpTC1N;
    TGraphErrors* gpTC2N;
    TGraphErrors* gpTC3N;
    TGraphErrors* gpTC4N;
    
    // Moments
    TGraphErrors* gpTM;   // Mean
    TGraphErrors* gpTVar; // Variance
    TGraphErrors* gpTS;   // Skewness
    TGraphErrors* gpTkap; // Kurtosis
    
    // Moments (Normalized to 0-5%)
    TGraphErrors* gpTMN;   // Mean
    TGraphErrors* gpTVarN; // Variance
    TGraphErrors* gpTSN;   // Skewness
    TGraphErrors* gpTkapN; // Kurtosis
    
    // -----**-----**-----**-----**-----**--N_{ch} Graphs--**-----**-----**-----**-----**
    
    // Cumulants (Absolute)
    TGraphErrors* gNchC1;
    TGraphErrors* gNchC2;
    TGraphErrors* gNchC3;
    TGraphErrors* gNchC4;
    
    // Cumulants (Normalized to 0-5%)
    TGraphErrors* gNchC1N;
    TGraphErrors* gNchC2N;
    TGraphErrors* gNchC3N;
    TGraphErrors* gNchC4N;
    
    // Moments
    TGraphErrors* gNchM;   // Mean
    TGraphErrors* gNchVar; // Variance
    TGraphErrors* gNchS;   // Skewness
    TGraphErrors* gNchkap; // Kurtosis
    
    // Moments (Normalized to 0-5%)
    TGraphErrors* gNchMN;   // Mean
    TGraphErrors* gNchVarN; // Variance
    TGraphErrors* gNchSN;   // Skewness
    TGraphErrors* gNchkapN; // Kurtosis
};

#endif // MOMENTHISTOGRAM_H
