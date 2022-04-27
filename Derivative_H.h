#ifndef CODE_H 
#define CODE_H

#include <algorithm>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <algorithm>
#include <array>
#include <iterator>  
#include <utility>

#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TList.h"
#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TF2.h"
#include "TH3F.h"
#include "TEnv.h"

using namespace std;

vector<float> MuBins = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42};
vector<float> PtBins = { 10.0, 30.0, 60.0, 120.0};
vector<float> EtaBins = {0.0, 0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.5};
vector<Int_t> NpvBins = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32};

vector<float> Config_EtaBins = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,3.0,3.075,3.15,3.25,3.35,3.45,3.6,3.8,4.1,4.5,4.9};
vector<float> Config_PtBins = { 10.0, 15.0, 20.0, 30.0, 60.0, 120.0, 200.0};
vector<Int_t> Config_NpvBins = {};
vector<float> Config_MuBins = {};

#endif