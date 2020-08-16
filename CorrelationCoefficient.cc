#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "TChain.h"

void CorrelationCoefficient()
{

    Float_t Mb_c, de;
    TChain *sig = new TChain("h1");
    sig->Add("code04singalMC.root");//input root file
    sig->SetBranchAddress("Mb_c", &Mb_c);
    sig->SetBranchAddress("de", &de);
    //sig->SetBranchAddress("etamass", &etamass);

    int total_n = (int)sig->GetEntries();

    double sumx = 0.;
    double sumy = 0.;
    double sumxy = 0.;
    double sumxsq = 0.;
    double sumysq = 0.;

    for (int k = 0; k < total_n; k++)
    {
        sig->GetEntry(k);

        sumx += mbc;
        sumy += de;
        sumxy += mbc * de;
        sumxsq += mbc * mbc;
        sumysq += de * de;
    }

    double numerator = total_n * sumxy - sumx * sumy;
    double denominator = sqrt((total_n * sumxsq - sumx * sumx) * (total_n * sumysq - sumy * sumy));

    cout << "Linear Correlation for all = " << (numerator / denominator) << endl;
    cout << "root Correlation for de and mbc = " << correlation(mbc, de) << endl;
}
