#include "TFile.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include <sstream>
#include <iostream>
#include <cstdio>
#include "TCanvas.h"
#include "TChain.h"
#include "THStack.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <cstdlib>

void TutorialCode()
{
    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");
    gROOT->ProcessLine("gPrintViaErrorHandler = kTRUE;");

    //Run 374810

    TChain *zdcchain = new TChain("zdcanalyzer/zdcdigi");
    TChain *hichain = new TChain("hiEvtAnalyzer/HiTree");
    TChain *hlttree = new TChain("hltanalysis/HltTree");
    TChain *skimchain = new TChain("skimanalysis/HltTree");
    TChain *trackchain = new TChain("PbPbTracks/trackTree");

    for (int x = 0 ; x < 1692 ; x++)
    {
        zdcchain->Add(Form("/store/data/Runs_2023/374810/HIPhysicsRawPrime01/HiForestMiniAOD_PbPb_HIPRP1_374810_%i.root/", x));
        hichain->Add(Form("/store/data/Runs_2023/374810/HIPhysicsRawPrime01/HiForestMiniAOD_PbPb_HIPRP1_374810_%i.root/", x));
        hlttree->Add(Form("/store/data/Runs_2023/374810/HIPhysicsRawPrime01/HiForestMiniAOD_PbPb_HIPRP1_374810_%i.root/", x));
        skimchain->Add(Form("/store/data/Runs_2023/374810/HIPhysicsRawPrime01/HiForestMiniAOD_PbPb_HIPRP1_374810_%i.root/", x));
        trackchain->Add(Form("/store/data/Runs_2023/374810/HIPhysicsRawPrime01/HiForestMiniAOD_PbPb_HIPRP1_374810_%i.root/", x));
    }
    zdcchain->AddFriend(hichain);
    zdcchain->AddFriend(hlttree);
    zdcchain->AddFriend(skimchain);
    zdcchain->AddFriend(trackchain);

    int zside[56];
    int section[56];
    int channel[56];

    float TS0[56];
    float TS1[56];
    float TS2[56];
    float TS3[56];
    float TS4[56];
    float TS5[56];

    for (int n = 0; n < 56; n++)
    {
        zside[n] = -1;
        section[n] = -1;
        channel[n] = -1;

        TS0[n] = -1.0;
        TS1[n] = -1.0;
        TS2[n] = -1.0;
        TS3[n] = -1.0;
        TS4[n] = -1.0;
        TS5[n] = -1.0;
    }

    float rpdsig[2][16] = {0.};

    for (int n = 0; n < 2; n++)
    {
        for (int m = 0; m < 16; m++)
        {
            rpdsig[n][m] = -1.0;
        }
    }
    
    float rpdx[4][2];
    float rpdy[4][2];

    for ( int n = 0; n < 4; n++)
    {
        for (int m = 0; m < 2; m++)
        {
            rpdx[n][m] = -5.;
            rpdy[n][m] = -5.;
        }
    }

    float rpdposx[2] = {0.,0.};
    float rpdposy[2] = {0.,0.};

    float rpdtotal[2] = {0.,0.};

    int nvtx = -1;
    std::vector<float> *zVtx = nullptr;
    int hibin = -1;
    double centrality = -1; 

    int MB = -1;

    int pprimaryVertexFilter = -1;
    int pclusterCompatibilityFilter = -1;
    int pphfCoincFilter2Th4 = -1;

    zdcchain->SetBranchAddress("zside",&zside);
    zdcchain->SetBranchAddress("section",&section);
    zdcchain->SetBranchAddress("channel",&channel);

    zdcchain->SetBranchAddress("chargefCTs0", &TS0);
    zdcchain->SetBranchAddress("chargefCTs1", &TS1);
    zdcchain->SetBranchAddress("chargefCTs2", &TS2);
    zdcchain->SetBranchAddress("chargefCTs3", &TS3);
    zdcchain->SetBranchAddress("chargefCTs4", &TS4);
    zdcchain->SetBranchAddress("chargefCTs5", &TS5);

    zdcchain->SetBranchAddress("nVtx",&nvtx);
    zdcchain->SetBranchAddress("hiBin",&hibin);
    zdcchain->SetBranchAddress("zVtx",&zVtx);

    //zdcchain->SetBranchAddress("HLT_HIMinimumBiasHF1ANDZDC1nOR_v1",&MB);

    zdcchain->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
    //zdcchain->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
    zdcchain->SetBranchAddress("pphfCoincFilter2Th4",&pphfCoincFilter2Th4);

    int NumEvents = zdcchain->GetEntries();
    cout << "The number of events in the root file is " << NumEvents << endl;
    NumEvents = 2000000;
    cout << "The number of events being analyzed in this root file is " << NumEvents << endl;
    for (int i = 0; i < NumEvents; i++)
    {
        zdcchain->GetEntry(i);

        if (i % 100000 == 0)
        {
            cout << "Event " << i << " is being processed" << endl;
        } 
             // Event selection
             if (fabs(zVtx->at(0)) > 15.0 || pphfCoincFilter2Th4 != 1 || pprimaryVertexFilter != 1)
             {
                continue;
             }

            for (int n = 0; n < 56; n++)
            {
                if (zside[n] == 1 && section[n] == 4)
                {
                    if (channel[n] == 1)  
                        { 
                            rpdsig[0][10] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 2)  
                        { 
                            rpdsig[0][6] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 3)  
                        { 
                            rpdsig[0][2] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 4)  
                        { 
                            rpdsig[0][14] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 5)  
                        { 
                            rpdsig[0][11] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 6)  
                        { 
                            rpdsig[0][7] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 7)  
                        { 
                            rpdsig[0][3] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 8)  
                        { 
                            rpdsig[0][15] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 9)  
                        { 
                            rpdsig[0][4] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 10) 
                        { 
                            rpdsig[0][8] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 11) 
                        { 
                            rpdsig[0][12] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 12) 
                        { 
                            rpdsig[0][0] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 13) 
                        { 
                            rpdsig[0][5] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 14) 
                        { 
                            rpdsig[0][9] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 15) 
                        { 
                            rpdsig[0][13] = (TS2[n] - TS1[n]);
                        }
                        if (channel[n] == 16) 
                        { 
                            rpdsig[0][1] = (TS2[n] - TS1[n]);
                        }
                    } // end RPD+

                    // RPD-

                    if (zside[n] == -1 && section[n] == 4)
                    {
                        if (channel[n] == 1)  
                            { 
                                rpdsig[1][10] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 2)  
                            { 
                                rpdsig[1][6] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 3)  
                            { 
                                rpdsig[1][2] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 4)  
                            { 
                                rpdsig[1][14] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 5)  
                            { 
                                rpdsig[1][11] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 6)  
                            { 
                                rpdsig[1][7] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 7)  
                            { 
                                rpdsig[1][3] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 8)  
                            { 
                                rpdsig[1][15] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 9)  
                            { 
                                rpdsig[1][4] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 10) 
                            { 
                                rpdsig[1][8] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 11) 
                            { 
                                rpdsig[1][12] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 12) 
                            { 
                                rpdsig[1][0] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 13) 
                            { 
                                rpdsig[1][5] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 14) 
                            { 
                                rpdsig[1][9] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 15) 
                            { 
                                rpdsig[1][13] = (TS2[n] - TS1[n]);
                            }
                            if (channel[n] == 16) 
                            { 
                                rpdsig[1][1] = (TS2[n] - TS1[n]);
                            }
                        } // end RPD-
                    }
        
            rpdtotal[0] = rpdsig[0][0] + rpdsig[0][1] + rpdsig[0][2] + rpdsig[0][3] + 
                            rpdsig[0][4] + rpdsig[0][5] + rpdsig[0][6] + rpdsig[0][7] + 
                            rpdsig[0][8] + rpdsig[0][9] + rpdsig[0][10] + rpdsig[0][11] + 
                            rpdsig[0][12] + rpdsig[0][13] + rpdsig[0][14] + rpdsig[0][15]; 

            rpdtotal[1] = rpdsig[1][0] + rpdsig[1][1] + rpdsig[1][2] + rpdsig[1][3] + 
                            rpdsig[1][4] + rpdsig[1][5] + rpdsig[1][6] + rpdsig[1][7] + 
                            rpdsig[1][8] + rpdsig[1][9] + rpdsig[1][10] + rpdsig[1][11] + 
                            rpdsig[1][12] + rpdsig[1][13] + rpdsig[1][14] + rpdsig[1][15]; 

            // =========================
            // Fill Histograms for signal distributions
            // =========================

        } // end events

        // Draw Histograms for signal distributions 
    } 