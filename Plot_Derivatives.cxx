#include "Derivative_H.h"
#include "new_try.cxx"

Double_t Get_Average_Slope(const vector<Double_t> &sl, const vector<Double_t> &we ){

    Double_t AV=0.0;
    Double_t WE=0.0;
    for(int t =0; t< sl.size();t++){
        WE+=we[t];
        AV+=we[t]*sl[t];
    }
    AV/=WE;
    return AV;
}

Double_t Get_RMS(const vector<Double_t> &sl, const vector<Double_t> &we , Double_t AV){

    Double_t RMS=0.0;
    Double_t WE=0.0;
    for(int t =0; t< sl.size();t++){
        WE+=we[t];
        
    }
    for(int t =0; t< sl.size();t++){
        Double_t w = we[t]/WE;
        RMS+=w*w*(sl[t]-AV)*(sl[t]-AV);
        
    }
    

    return sqrt(RMS);
}
Double_t Get_Averaged_Error(const vector<Double_t> &er, const vector<Double_t> &we ){

    Double_t ER=0.0;
    Double_t WE=0.0;
    for(int t =0; t< er.size();t++){
        WE+=we[t];        
    }
    for(int t =0; t< er.size();t++){
        Double_t w = we[t]/WE;
        ER+=w*w*er[t]*er[t];
    }
    

    return sqrt(ER);
}
Double_t Get_MAX_Error(const vector<Double_t> &er){

    Double_t ER=0.0;
    
    
    for(int t =0; t< er.size();t++){
        if(abs(er[t])>ER) ER = abs(er[t]);
    }
    

    return ER;
}
void Plot_Derivatives(){
    TProfile * TProfile_Pt_vs_mu_binned [ (EtaBins.size() -1) * (PtBins.size()) *( NpvBins.size()-1 ) ];
    
    
    //      /afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/saved_TProfiles.root  NO corr
    //      /afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/project_derivative/saved_TProfiles_area_coor.root  AREA correction
    TFile * loaded_TProfiles = new TFile("/afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/project_derivative/saved_TProfiles_area_coor.root","read");
    for ( int i=0; i< EtaBins.size()-1; i++ ){
        for ( int j=0; j< PtBins.size(); j++ ){
            for ( int k=0; k< NpvBins.size()-1; k++ ){
                int iter = i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k;
                string name_str = Create_Name_eta_pt_NPV_mu("TProfile_Pt_VS_Mu_Eta_",i, j, k, 0, iter );
                TProfile_Pt_vs_mu_binned[iter] = (TProfile*)loaded_TProfiles->Get(From_String_To_Char_Array(name_str));
                
            }
        }
    }
    TProfile_Pt_vs_mu_binned[0]->Print();
    

    TCanvas* Averaged_pt_canvas[ (EtaBins.size() -1) * (PtBins.size()) *( NpvBins.size()-1 )];
    TF1 * averaged_linear_fit[ (EtaBins.size() -1) * (PtBins.size()) *( NpvBins.size()-1 ) ];
    for ( int i=0; i< EtaBins.size()-1; i++ ){
        for ( int j=0; j< PtBins.size(); j++ ){
            for ( int k=0; k< NpvBins.size()-1; k++ ){
                int iter = i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k;
                vector<int> bounaries = Boundaries_Of_Fit((*TProfile_Pt_vs_mu_binned[iter]));
                
                string name_str = Create_Name_eta_pt_NPV_mu("averaged_canvas",i, j, k, 0, iter );
                Averaged_pt_canvas[iter] = new TCanvas(From_String_To_Char_Array(name_str), From_String_To_Char_Array(name_str), 800, 800);
                Averaged_pt_canvas[iter ]->cd();
                
                TProfile_Pt_vs_mu_binned[iter]->SetLineColor(kBlue);
                TProfile_Pt_vs_mu_binned[iter]->Draw("E1");
                string name_fit_name_str =    Create_Name_eta_pt_NPV_mu("averaged_linear_fit_",i, j, k, 0, iter );
                //cout<<From_String_To_Char_Array(name_fit_name_str)<<endl;
                
                averaged_linear_fit[iter ] = new TF1(From_String_To_Char_Array(name_fit_name_str) ,"[0]+x*[1]",MuBins[0],MuBins[MuBins.size()-1]);
                averaged_linear_fit[iter ]->SetParameter(0,0.0);
                averaged_linear_fit[iter ]->SetParameter(1,0.0);
                
                //TProfile_Pt_vs_mu_binned[(i-1)*PtBins.size() + j]->Fit(averaged_linear_fit[(i-1)], "F", "L",MuBins[0], 45.0); //
                if(Total_TProfile_Entries(*(TProfile_Pt_vs_mu_binned[iter])) != 0){
                    if(bounaries[0]>=0 && bounaries[1]>0 && bounaries[1] - bounaries[0] - 1 >=2){
                        TProfile_Pt_vs_mu_binned[iter]->Fit(averaged_linear_fit[iter], "M","A",MuBins[bounaries[0]+1],MuBins[bounaries[1]-1]); 
                    }
                        
                    else{
                        averaged_linear_fit[iter ]->SetParameter(0,0.0);
                        averaged_linear_fit[iter ]->SetParameter(1,0.0);
                    }
                }
                
                cout<<Total_TProfile_Entries((*TProfile_Pt_vs_mu_binned[iter]))<<"\n";
                cout<<iter <<endl;
                
                //cout<< EtaBins[i] << ": "<<  averaged_linear_fit[(i-1)*PtBins.size() + j]->GetParameter(1) <<endl;
                averaged_linear_fit[iter]->SetLineColor(kRed);
                averaged_linear_fit[iter]->SetLineWidth(3);
                Averaged_pt_canvas[iter]->cd();
                averaged_linear_fit[iter]->Draw("SAME");
                
            
                name_str = Create_Name_eta_pt_NPV_mu("/afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/new_canvases/averaged_",i, j, k, 0, iter ,"_pt_vs_mu.png" );
                //Averaged_pt_canvas[iter]->SaveAs(From_String_To_Char_Array(name_str));
                cout<<"\n";
                /*
                */
            }
        }
    }

    
    TH1D * THistos[PtBins.size()];

    for ( int j=0; j< PtBins.size(); j++ ){
        string name = Create_Name_eta_pt_NPV_mu("th1d",0, j, 0, 0, j );
        THistos[j] = new TH1D(From_String_To_Char_Array(name),From_String_To_Char_Array(name),EtaBins.size()-1, &EtaBins[0]);
    }
    
    for ( int i=0; i< EtaBins.size()-1; i++ ){

        for ( int j=0; j< PtBins.size(); j++ ){
            vector <Double_t> slopes={};
            vector <Double_t> weights={};
            vector <Double_t> errors={};
            for ( int k=0; k< NpvBins.size()-1; k++ ){
                int iter = i + ( EtaBins.size() - 1 )*j +  ( EtaBins.size() - 1 )*( PtBins.size() )*k;
                
                slopes.push_back(averaged_linear_fit[iter]->GetParameter(1));
                weights.push_back(Total_TProfile_Entries(*(TProfile_Pt_vs_mu_binned[iter])));
                errors.push_back(averaged_linear_fit[iter]->GetParError(1));
            }
            cout<<"ETABIN: "<<i<<" PTBIN: "<<j<<"\n";
            for(int t =0; t< slopes.size();t++ ){
                cout<< slopes[t] << "; ";
            }
            cout<<"\n";
            cout<<"Average: "<< Get_Average_Slope(slopes, weights)<<"\n";
            cout<<"RMS: "<<Get_RMS(slopes,weights, Get_Average_Slope(slopes, weights))<<"\n";
            cout<<"Error: "<<Get_Averaged_Error(errors,weights)<<"\n";
            cout<<"MAX Error: "<<Get_MAX_Error(errors)<<"\n";
            for(int t =0; t< slopes.size();t++ ){
                cout<< errors[t] << "; ";
            }
            cout<<"\n";
            cout<<"\n";
            THistos[j]->SetBinContent(i+1, Get_Average_Slope(slopes, weights));
            THistos[j]->SetBinError(i+1,Get_Averaged_Error(errors,weights));
            cout<< THistos[j]->GetBinWidth(i+1)<<"\n";
        }

    }
   /**/
    TCanvas* TH1d_pt_canvas[ PtBins.size()];
    for ( int j=0; j< PtBins.size(); j++ ){
        string name = Create_Name_eta_pt_NPV_mu("canvas",0, j, 0, 0, j );
        TH1d_pt_canvas[j] = new TCanvas(From_String_To_Char_Array(name), From_String_To_Char_Array(name), 800, 800);
        TH1d_pt_canvas[j ]->cd();
        THistos[j]->Draw();
        THistos[j]->Print();
        name = Create_Name_eta_pt_NPV_mu("/afs/cern.ch/work/d/dtimoshy/RC/Draw_Utensils/new_canvases/th1d_",0, j, 0, 0, j ,"_dpt_vs_eta.png" );
        TH1d_pt_canvas[j]->SaveAs(From_String_To_Char_Array(name));
    }
     
    for ( int i=0; i< EtaBins.size(); i++ ){
        cout<< EtaBins[i]<<";";
        //cout<< &EtaBins[0];
    }
    cout<<"\n";
    
    loaded_TProfiles->Close();
    
    
    return;
}