#include "Derivative_H.h"



bool BOOL_exists_file (const string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}
char * From_String_To_Char_Array( string & name){
    char * char_name[500];
    for (int i =0; i < name.size();i++){
        //cout<<name.at(i);
        char_name[i] = & name.at(i);
    }
    
    return (*char_name);
}

// This function gives the bin number in which the event by value of eta should be stored. 
// Can be used for general stuff

int Get_Bin_Of_Absolute_Value( vector<float> & etabin , float &eta ){
    for (int i = 0 ; i< etabin.size(); i++ ){
        if (abs(eta) <= etabin[i]){
            if (i == 0)
                return -1;
            else 
                return i;
        }
    }
    return -1;
}
template <typename T> 
int Get_Bin_Of_Absolute_Value_T( vector<T> & etabin ,  T &eta ){
    for (int i = 0 ; i< etabin.size(); i++ ){
        if (abs(eta) <= etabin[i]){
            if (i == 0)
                return -1;
            else 
                return i;
        }
    }
    return -1;
}

int Get_Bin_Of_Value( vector<float> & etabin , float &eta ){
    for (int i = 0 ; i< etabin.size(); i++ ){
        if ( eta <= etabin[i]){
            if (i == 0)
                return -1;
            else 
                return i;
        }
    }
    return -1;
}

// Returns boundaries of filled elements in TProfile
Double_t Total_TProfile_Entries(TProfile & TProf){
    Double_t total = 0.0;
    for (int i =1 ; i < TProf.GetNbinsX()+1; i++) {
        total += TProf.GetBinEntries(i);
        
        
    }
    return total;
}

vector<int> Boundaries_Of_Fit(TProfile & TProf){
    vector<int> a;
    cout<< TProf.GetNbinsX()<<endl;
    cout<< TProf.GetNbinsY()<<endl;
    int x_min=0;
    int x_max=0;
    for (int i =1 ; i < TProf.GetNbinsX()+1; i++) {
        if ( TProf.GetBinContent( i )>0){
            if(x_min==0)
                x_min=i;
            x_max=i;
        }
        cout<< TProf.GetBinContent(i)<<"; ";
        
    }
    
    cout<<"\n";
    a.push_back(x_min);
    a.push_back(x_max);
    cout<<"First: "<<a[0]<<" Second: "<<a[1]<<"\n";
    return a;

}
template <typename T , typename Y> 
int Get_Bin_Of_Absolute_Value_T( vector<T> & etabin ,  Y &eta ){
    for (int i = 0 ; i< etabin.size(); i++ ){
        if (abs(eta) <= etabin[i]){
            if (i == 0)
                return -1;
            else 
                return i;
        }
    }
    return -1;
}
template <typename T , typename Y,  typename U, typename I, typename O> 
string Create_Name_eta_pt_NPV_mu(string custom, T eta, Y pt, U NPV, I mu, O n, string end = "" ){
    string out = "";
    out = custom;
    out += "_eta_"+to_string((int)(eta));
    out += "_pt_"+to_string((int)pt);
    out += "_NPV_"+to_string((int)NPV);
    out += "_mu_"+to_string((int)mu);
    out += "_n_"+to_string((int)n);
    out += end;
    return out;
}

Double_t DeltaPtResidual(Double_t pt, Double_t z, Double_t f,Double_t s ){
    return z + f*pt+s*TMath::Log(pt);
}
Double_t Get_Correction(TList * CoeffLists, Double_t pt, int eta_bin, int NPV_bin, int mu_bin){
    string zero_order_name = "intersections";
    string first_order_name = "slopes";
    string second_order_name = "SecondOrderCoeff";
    //CoeffLists->Print();
    TList * intercepts = (TList *) CoeffLists->FindObject(From_String_To_Char_Array(zero_order_name));	
    
	TList * slopes = (TList *) CoeffLists->FindObject(From_String_To_Char_Array(first_order_name));	
    
	TList * secondOs = (TList *) CoeffLists->FindObject(From_String_To_Char_Array(second_order_name));	

    Double_t coeff_i = 0.0 ;
	Double_t coeff_s = 0.0 ;
	Double_t coeff_SO = 0.0 ;

    string name_of_th2d = "intersections_eta"+to_string(eta_bin);
	TH2D * TH2D_intercept = (TH2D*)intercepts->FindObject(From_String_To_Char_Array(name_of_th2d));
	name_of_th2d = "slopes_eta"+to_string(eta_bin);
    TH2D * TH2D_slope = (TH2D*)slopes->FindObject(From_String_To_Char_Array(name_of_th2d));
    name_of_th2d = "intersections_eta"+to_string(eta_bin);
    TH2D * TH2D_secondO = (TH2D*)secondOs->FindObject(From_String_To_Char_Array(name_of_th2d));

    coeff_i = TH2D_intercept->GetBinContent(NPV_bin,mu_bin);
    coeff_s = TH2D_slope->GetBinContent(NPV_bin,mu_bin);
    coeff_SO = TH2D_secondO->GetBinContent(NPV_bin,mu_bin);
    
    delete TH2D_intercept, TH2D_slope, TH2D_secondO;
    
    return DeltaPtResidual(pt,coeff_i,coeff_s,coeff_SO );
}

vector<double> ComputeOffsets(vector<double> etaBins, vector<double> term){
  std::vector <double> offset;
  offset.push_back(term.at(0));  
  for(unsigned int i=1; i<etaBins.size();++i){
    double espace = etaBins.at(i)-etaBins.at(i-1);
    double ioffs = offset.at(i-1)+term.at(i)*espace;
    offset.push_back(ioffs);
  }
  return offset;
}

double interpolation (vector<double> etaBins, vector<double> term, vector<double> offset, Float_t eta){
  double correction=0;
  //  std::cout<<"eta : "<<eta<<std::endl;
  for(unsigned int i=0; i<(etaBins.size()-1);++i){
    if (eta >= etaBins.at(i) && eta < etaBins.at(i+1)){
      correction = offset.at(i) + (eta-etaBins.at(i))*term.at(i+1); // computing value of the linear function at bin i for a given value of eta
      //  std::cout<<"inside interpolation if"<<std::endl;
      break;
    }
  }
  return correction;      
}


void Create_TProfiles(){
    
    cout<<"Welcome Aristocrat!"<<endl;

    



    
    string str;
    ifstream infile ;
    vector <string> paths_to_ttrees ;
    infile.open("/afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/List_of_ttrees.txt");

    while(!infile.eof()){
        getline(infile,str);
        paths_to_ttrees.push_back(str);
        if(BOOL_exists_file(str))
            cout<<"FILE EXISTS: "<<str<<"\n\n";
        else
            cout<<"ERROR FILE DOES NOT EXIST: "<<str<<"\n\n";
    }
    infile.close();

    

    
    
    vector<Float_t>* pt; //true pt    
    Int_t           NPV;
    Float_t         mu;
    Double_t        weight_tot;
    Float_t        weight;
    std::vector<float>* jet_eta; //jet eta
    vector<Float_t>* pt_true;
    vector<Float_t>* jet_area;
    Float_t         rho;

    TBranch* b_pt;
    TBranch* b_eta;
    TBranch* b_NPV;
    TBranch* b_mu;
    TBranch* b_weight;
    TBranch* b_weight_tot;
    TBranch* b_pt_true;
    TBranch* b_jet_area;
    TBranch* b_rho;
    
    
    //vector<float> PtBins = { 0.0, 200.0};
    //vector<float> EtaBins = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,3.0,3.075,3.15,3.25,3.35,3.45,3.6,3.8,4.1,4.5,4.9};
    //vector<float> EtaBins = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0,1.2, 1.4,1.6,1.8,2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 3.9, 4.2, 4.5};
    //vector<Int_t> NpvBins = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74};
    //vector<float> PtBins = { 10.0, 15.0, 20.0, 30.0, 60.0, 120.0, 200.0};


    
   
    

    
    
    // TProfile in which values are stored
    // They are organized, such that each event of bin in eta and in pt is saved in corresponding 
    // number of TProfile array. 
    // If event has eta bin I and pt bin J, than it is stored in: (I-1)*PtBins.size() + J
    vector<int> check;
    TProfile * TProfile_Pt_vs_mu_binned [ (EtaBins.size() -1) * (PtBins.size()) *( NpvBins.size()-1 ) ];
    for ( int i=0; i< EtaBins.size()-1; i++ ){
        for ( int j=0; j< PtBins.size(); j++ ){
            for ( int k=0; k< NpvBins.size()-1; k++ ){
                int iter = i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k ;
                string name_str = Create_Name_eta_pt_NPV_mu("TGraph_Pt_VS_Mu",i, j, k, 0, iter );
                
                check.push_back(iter);
                
                name_str = Create_Name_eta_pt_NPV_mu("TProfile_Pt_VS_Mu_Eta_",i, j, k, 0, iter,  ";#mu;p_{T}");
                
                TProfile_Pt_vs_mu_binned[ i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k ] = new TProfile(From_String_To_Char_Array(name_str),From_String_To_Char_Array(name_str), MuBins.size(), MuBins[0],MuBins[MuBins.size()-1],  0.0,10000.0 );
                TProfile_Pt_vs_mu_binned[ i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k ]->Sumw2();
                cout<< From_String_To_Char_Array(name_str)<<endl;
            }
        }
    }
    
    
    
    

    string path_to_parameters = "/afs/cern.ch/work/d/dtimoshy/RC/saved/residualCalibParameters.root";	
	TFile *f = new TFile(From_String_To_Char_Array(path_to_parameters));
	TList *CoeffLists = (TList*)f->Get("param3D");
    
    
    for (int i = 0 ; i <=  paths_to_ttrees.size()-1 ; i++){
        TFile* FILE_TO_TTREE = new TFile(From_String_To_Char_Array(paths_to_ttrees[i]),"read");
        cout<< "\n\n Processing file: "<<paths_to_ttrees[i]<<"\n";
        TTree* Tree = (TTree*) FILE_TO_TTREE->Get("IsolatedJet_tree");
        Tree->SetBranchAddress("jet_ConstitPt", &pt, &b_pt);
        Tree->SetBranchAddress("weight_tot", &weight_tot, &b_weight_tot);
        Tree->SetBranchAddress("weight", &weight, &b_weight);
        Tree->SetBranchAddress("actualInteractionsPerCrossing", &mu, &b_mu);
        Tree->SetBranchAddress("jet_eta", &jet_eta, &b_eta);
        Tree->SetBranchAddress("jet_true_pt", &pt_true, &b_pt_true);
        Tree->SetBranchAddress("jet_ActiveArea4vec_pt", &jet_area, &b_jet_area);
        Tree->SetBranchAddress("rho", &rho, &b_rho);
        Tree->SetBranchAddress("NPV", &NPV, &b_NPV);
        //recojet_pt->at(j)-pt->at(j)-jet_area->at(j)*rho*0.001
        for(int ientry=0;ientry<Tree->GetEntries()/10000;ientry++){ 
            Tree->GetEntry(ientry);
            for (int jet_iter = 0; jet_iter < pt->size(); jet_iter++  ){
                //cout<< pt->at(jet_iter)<<" ";
                
                int eta_bin = Get_Bin_Of_Absolute_Value( EtaBins, jet_eta->at(jet_iter) );
                int pt_bin = Get_Bin_Of_Absolute_Value( PtBins, pt_true->at(jet_iter) );
                int NPV_bin = Get_Bin_Of_Absolute_Value_T( NpvBins, NPV);
                //cout<< eta_bin << " "<< pt_bin<<endl;
                //cout<< jet_eta->at(jet_iter)  << " " <<  pt_true->at(jet_iter)<<endl ;
                
                int iter = ( eta_bin -1 ) + ( EtaBins.size() - 1 )*( pt_bin -1 ) +  ( EtaBins.size() - 1 )*( PtBins.size() )*( NPV_bin -1 );
                int iter_inclusive = ( eta_bin -1 ) + ( EtaBins.size() - 1 )*( PtBins.size() -1 ) +  ( EtaBins.size() - 1 )*( PtBins.size() )*( NPV_bin -1 );
                if (eta_bin > 0 && pt_bin > 0 && NPV_bin > 0){
                    //Double_t pt_area = pt->at(jet_iter) - jet_area->at(jet_iter)*rho;
                    //Double_t DeltaResidual = Get_Correction(CoeffLists, pt_area,  eta_bin-1,  npv_bin-1,  mu_bin-1);

                    //if (DeltaResidual == DeltaResidual && pt_area - DeltaResidual > 0){
                        //TProfile_Pt_vs_mu_binned [iter]->Fill( mu, pt_area - DeltaResidual , weight  );
                    //}
                    TProfile_Pt_vs_mu_binned [iter]->Fill( mu, pt->at(jet_iter) , weight  );
                    //TProfile_Pt_vs_mu_binned [iter]->Fill( mu, pt->at(jet_iter) - jet_area->at(jet_iter)*rho, weight  );
                    //TProfile_Pt_vs_mu_binned [iter]->Fill( mu, pt_true->at(jet_iter) , weight  );
                }
                if (eta_bin > 0 && NPV_bin > 0){

                    //Double_t pt_area = pt->at(jet_iter) - jet_area->at(jet_iter)*rho;
                    //Double_t DeltaResidual = Get_Correction(CoeffLists, pt_area,  eta_bin-1,  npv_bin-1,  mu_bin-1);
                    //if (DeltaResidual == DeltaResidual && pt_area - DeltaResidual > 0){
                        //TProfile_Pt_vs_mu_binned [iter]->Fill( mu, pt_area - DeltaResidual , weight  );
                    //}

                    TProfile_Pt_vs_mu_binned [iter_inclusive]->Fill( mu, pt->at(jet_iter) , weight  );
                    //TProfile_Pt_vs_mu_binned [iter_inclusive]->Fill( mu, pt->at(jet_iter) - jet_area->at(jet_iter)*rho , weight  );
                    //TProfile_Pt_vs_mu_binned [iter]->Fill( mu, pt_true->at(jet_iter) , weight  );
                }
                if ( eta_bin == 1 && pt_bin == 1){
                    //cout<<( eta_bin - 1 ) * PtBins.size() + pt_bin-1<<endl;
                    //cout<< eta_bin << " "<< pt_bin<<endl;
                    //cout<< jet_eta->at(jet_iter)  << " " <<  pt_true->at(jet_iter)<<endl ;
                }
                    
                


            }
          
        }
        cout<<Tree->GetEntries()<<"\n";
        FILE_TO_TTREE->Close();
    }
    f->Close();
    string name_to_save = "/afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/project_derivative/saved_TProfiles_NO_corr.root";
    TFile * saved_TProfiles = new TFile(From_String_To_Char_Array(name_to_save),"recreate");
    for ( int i=0; i< EtaBins.size()-1; i++ ){
        for ( int j=0; j< PtBins.size(); j++ ){
            for ( int k=0; k< NpvBins.size()-1; k++ ){
                int iter = i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k ;
                string name_str = Create_Name_eta_pt_NPV_mu("TProfile_Pt_VS_Mu_Eta_",i, j, k, 0, iter );
                TProfile_Pt_vs_mu_binned[iter]->Write(From_String_To_Char_Array(name_str),TObject::kSingleKey);
                
            }
        }
    }
    saved_TProfiles->Close();


    
}



