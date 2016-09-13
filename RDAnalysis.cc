//RDAnalysis.C
//****************************************************
//This program computes the RDParameter for the CMS 
//run numbers listed in certainruns[]
//****************************************************

#include "MyParser.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TTree.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TDatime.h"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

// runs with pots insertion that you want to analyze
unsigned int certainruns[] = {273725,273728,274094,274172,274198,274199,274200,274241,274244,274284,274286,274387,274388,274422,274440,274441,
  274442,274443,274955,274958,274968,274969,274970,274969,274999,275000,275001,275066,275074,275124,275125,275282,275283,275290,275292,275293,275309,
  275310,275311,275319,275338,275344,275345,275370,275371,275375,275376,275768,275772,275774,275776,275777,275778,275782,275783,275828,275829,
  275832,275833,275834,275836,275837,275847,275890,275911,275912,275913,275918,275920,
279766,
279794,
279823,
279841,
279849,
279931,
279966,
279975,
280018,
280191,
280330,
280385
}; 

// runs with data, without pots insertion (read from fillreport)
unsigned int excludedruns[] = {274100,274102,274103,274104,274105,274106,274107,274108,274142,274146,274157,274159,274160,274161,274250,
  274251,274314,274315,274316,274317,274318,274319,274335,274336,274337,274338,274339,274344,274345,274382,274966,275326};

const int ncertainruns = sizeof(certainruns)/sizeof(unsigned int);
const int nexcludedruns = sizeof(excludedruns)/sizeof(unsigned int);

vector <double> luminosities_int(ncertainruns,0);
vector <double> ycertainruns_nr_rp(ncertainruns,0);
vector <double> ycertainruns_fr_rp(ncertainruns,0);
vector <double> eycertainruns_nr_rp(ncertainruns,0);
vector <double> eycertainruns_fr_rp(ncertainruns,0);
vector <double> temp_sum_normalized_nr_rp(ncertainruns,0);
vector <double> temp_sum_normalized_fr_rp(ncertainruns,0);
vector <double> interaction_BX_nr_rp(ncertainruns,0);
vector <double> interaction_BX_fr_rp(ncertainruns,0);
vector <double> start_times(ncertainruns,0);
// vector <double> stop_times(ncertainruns,0);
vector <TDatime> datimes(ncertainruns,TDatime());

//----------------------------------------------------------------------------------------------------

int main()
{
	TFile *f_out = new TFile("output.root", "recreate");
  
	ifstream input_file("runinfo.csv");
	vector<string> lines;
	vector<vector<string>> fields;
	LineParser(input_file, lines);
	CSVLinesParser(lines, fields);

	sort(certainruns, certainruns+ncertainruns);
	sort(excludedruns, excludedruns+nexcludedruns);

	for (int i = 0; i < ncertainruns; ++i)
	{
		double temp_sum_luminosities_int=0;
		for(vector<vector<string>>::reverse_iterator lit = fields.rbegin(); lit != fields.rend(); ++lit)
		{
			if (lit->size() < 8)
			{
				printf("ERROR: less than 8 fields on the line. Skipping.\n");
				continue;
			}

			unsigned int f_run = atoi(lit->at(0).c_str());
			double f_lumi = atof(lit->at(2).c_str()) / 1000000;

			TDatime datetime(atoi(lit->at(7).substr(0,4).c_str()), atoi(lit->at(7).substr(5,2).c_str()),
				atoi(lit->at(7).substr(8,2).c_str()), atoi(lit->at(7).substr(11,2).c_str()),
				atoi(lit->at(7).substr(14,2).c_str()), atoi(lit->at(7).substr(17,2).c_str()));

			bool runToExclude = false;
			for (unsigned int er : excludedruns)
			{
				if (f_run == er)
				{
					runToExclude = true;
					break;
				}
			}

			if (!runToExclude)
				temp_sum_luminosities_int += f_lumi;

			if (f_run == certainruns[i])
			{
				luminosities_int[i] = temp_sum_luminosities_int;
				datimes[i] = datetime;
				start_times[i] = datimes[i].Convert();
			}
		}
	}
  
  for (int i=0;i<ncertainruns;i++){
    stringstream ss;
    ss << "Runs/DQM_V0001_CTPPS_R000" << certainruns[i] << ".root"; //name of the root file
//     cout << ss.str() << endl;
    TFile *file = TFile::Open((ss.str()).c_str());
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/TrackingStrip/sector 45/station 210/nr_hr/track XY profile"; //path to the histogram
    for (int binx=12;binx<19;binx++){ //compute integral over the selected area
      for (int biny=50;biny<52;biny++){
	if (certainruns[i]<274244){
	ycertainruns_nr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
	}
	else{
	ycertainruns_nr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny);
	}
      }
    }
    
    double temp_sum = 0;
    
    for (int binx=24;binx<37;binx++){ //compute integral over the reference area
      for (int biny=33;biny<38;biny++){
	if (certainruns[i]<274244){
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
	}
	else{
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny); 
	}
      }
    }    
//     ycertainruns_nr_rp.at(i) /= (((TH2F*) file->Get(ss.str().c_str()))->GetEntries() - ycertainruns_nr_rp.at(i));
    eycertainruns_nr_rp.at(i) = ycertainruns_nr_rp.at(i) * sqrt(1./ycertainruns_nr_rp.at(i) + 1./temp_sum) / temp_sum;
    ycertainruns_nr_rp.at(i) /= temp_sum;
//     eycertainruns_nr_rp.at(i) = sqrt(ycertainruns_nr_rp.at(i))/luminosities[i];
//     ycertainruns_nr_rp.at(i) /= luminosities[i];
//     cout << ycertainruns_nr_rp.at(i) << endl;
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/events per BX";
    temp_sum_normalized_nr_rp.at(i) = temp_sum;
    interaction_BX_nr_rp.at(i) = ((TH2F*) file->Get(ss.str().c_str()))->GetEntries();
//     etemp_sum_normalized_nr_rp.at(i) = temp_sum_normalized_nr_rp.at(i)*sqrt(1./((TH2F*) file->Get(ss.str().c_str()))->GetEntries() + 1./temp_sum);
    ss.str(""); //resets string
}

  
//   calculates parameter for the 45 fr rp with all the "good" runs
// cout << ncertainruns << endl;
  for (int i=0;i<ncertainruns;i++){
    stringstream ss;
    ss << "Runs/DQM_V0001_CTPPS_R000" << certainruns[i] << ".root"; //name of the root file
//     cout << ss.str() << endl;
    TFile *file = TFile::Open((ss.str()).c_str());
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/TrackingStrip/sector 45/station 210/fr_hr/track XY profile"; //path to the histogram
    for (int binx=18;binx<25;binx++){ //compute integral over the selected area
      for (int biny=40;biny<42;biny++){
	if (certainruns[i]<274244){
	ycertainruns_fr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
	}
	else{
	ycertainruns_fr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny); 
	}
      }
    }
    
    double temp_sum = 0;
    for (int binx=30;binx<43;binx++){ //compute integral over the reference area
      for (int biny=23;biny<28;biny++){
	if (certainruns[i]<274244){
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
	}
	else{
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny); 
	}
      }
    }
//     ycertainruns_fr_rp.at(i) /= (((TH2F*) file->Get(ss.str().c_str()))->GetEntries() - ycertainruns_fr_rp.at(i));
    eycertainruns_fr_rp.at(i) = ycertainruns_fr_rp.at(i) * sqrt(1./ycertainruns_fr_rp.at(i) + 1./temp_sum) / temp_sum;
    ycertainruns_fr_rp.at(i) /= temp_sum;
//     eycertainruns_fr_rp.at(i) = sqrt(ycertainruns_fr_rp.at(i))/luminosities[i];
//     ycertainruns_fr_rp.at(i) /= luminosities[i];   
//     cout << ycertainruns_nr_rp.at(i) << endl;
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/events per BX";
    temp_sum_normalized_fr_rp.at(i) = temp_sum;
    interaction_BX_fr_rp.at(i) = ((TH2F*) file->Get(ss.str().c_str()))->GetEntries();
//     etemp_sum_normalized_fr_rp.at(i) = temp_sum_normalized_fr_rp.at(i)*sqrt(1./((TH2F*) file->Get(ss.str().c_str()))->GetEntries() + 1./temp_sum);
    ss.str(""); //resets string
  }

// graphics part
//c1,c2 plots of parameter vs total integrated luminosity
//check_ref_area plots of correlation between hits in the reference area and # interactions per BX
//c1_time,c2_time plots of parameter vs time

	gDirectory = f_out;

  TCanvas *c1 = new TCanvas("c1","RDAnalysis_rp_45_210_nr_hr",700,500);
  c1->SetGrid();
  c1->SetLogy();
//   TGraphErrors *gr1 = new TGraphErrors(ncertainruns,certainruns,&(ycertainruns_nr_rp[0]),NULL,&(eycertainruns_nr_rp[0]));
  TGraphErrors *gr1 = new TGraphErrors(ncertainruns,&(luminosities_int[0]),&(ycertainruns_nr_rp[0]),NULL,&(eycertainruns_nr_rp[0]));
  TF1 *l1 = new TF1("l1","51.9391304",0,7);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(0.6);
  gr1->SetTitle("RDAnalysis_rp_45_210_nr_hr");
//   gr1->GetXaxis()->SetTitle("Run");
  gr1->GetXaxis()->SetTitle("Integrated Luminosity (fb^{-1})");  
  gr1->GetYaxis()->SetTitle("Parameter"); 
  gr1->Draw("AP");
  l1->Draw("same");

  c1->Write("45_210_nr_hr, lumi");
  
  TCanvas *c2 = new TCanvas("c2","RDAnalysis_rp_45_210_fr_hr",700,500);
  c2->SetGrid();
  c2->SetLogy();
//   TGraphErrors *gr2 = new TGraphErrors(ncertainruns,certainruns,&(ycertainruns_fr_rp[0]),NULL,&(eycertainruns_fr_rp[0])); 
  TGraphErrors *gr2 = new TGraphErrors(ncertainruns,&(luminosities_int[0]),&(ycertainruns_fr_rp[0]),NULL,&(eycertainruns_fr_rp[0]));   
  TF1 *l2 = new TF1("l2","55.1898",0,7);  
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(0.6);
  gr2->SetTitle("RDAnalysis_rp_45_210_fr_hr");
//   gr2->GetXaxis()->SetTitle("Run");
  gr2->GetXaxis()->SetTitle("Integrated Luminosity (fb^{-1})");  
  gr2->GetYaxis()->SetTitle("Parameter");
  gr2->SetMaximum(65.);
  gr2->Draw("AP");
  l2->Draw("same");
  
  c2->Write("45_210_fr_hr, lumi");
  
//Graphics to check if the ratio between reference area hits and events per bunch crossing is constant
  TCanvas *check_ref_area = new TCanvas("check_ref_area","rp_45_210");
  check_ref_area->Divide(2);
  TPad* check_ref_area_1 = (TPad*)(check_ref_area->GetPrimitive("check_ref_area_1"));
  TPad* check_ref_area_2 = (TPad*)(check_ref_area->GetPrimitive("check_ref_area_2"));
  check_ref_area_1->cd();
  TGraph *check_ref_area_fr = new TGraph(ncertainruns,&(temp_sum_normalized_fr_rp[0]),&(interaction_BX_fr_rp[0]));   
  check_ref_area_fr->SetMarkerStyle(21);
  check_ref_area_fr->SetMarkerSize(0.6);
  check_ref_area_fr->SetTitle("Reference area vs interactions correlation 45_210_fr_hr");
  check_ref_area_fr->GetXaxis()->SetTitle("Hits in the reference area");
  check_ref_area_fr->GetYaxis()->SetTitle("Total number of triggers");
  check_ref_area_fr->GetYaxis()->SetTitleOffset(1.6);
  check_ref_area_fr->Draw("AP");
  check_ref_area_2->cd();
  TGraph *check_ref_area_nr = new TGraph(ncertainruns,&(temp_sum_normalized_nr_rp[0]),&(interaction_BX_nr_rp[0]));   
  check_ref_area_nr->SetMarkerStyle(21);
  check_ref_area_nr->SetMarkerSize(0.6);
  check_ref_area_nr->SetTitle("Reference area vs interactions correlation 45_210_nr_hr");
  check_ref_area_nr->GetXaxis()->SetTitle("Hits in the reference area");
  check_ref_area_nr->GetYaxis()->SetTitle("Total number of triggers");
  check_ref_area_nr->GetYaxis()->SetTitleOffset(1.6);
  check_ref_area_nr->Draw("AP");

  check_ref_area->Write("45_210, check_ref_area");
  
  TCanvas *c1_time = new TCanvas("c1_time","RDAnalysis_rp_45_210_nr_hr_time",700,500);
  c1_time->SetGrid();
  c1_time->SetLogy();
  TGraphErrors *gr1_time_start = new TGraphErrors(ncertainruns,&(start_times[0]),&(ycertainruns_nr_rp[0]),NULL,&(eycertainruns_nr_rp[0]));
  TF1 *l1_time = new TF1("ref","51.9391304",start_times.front(),start_times.back());
  gr1_time_start->SetMarkerStyle(21);
  gr1_time_start->SetMarkerSize(0.6);
  gr1_time_start->SetTitle("RDAnalysis_rp_45_210_nr_hr");
  gr1_time_start->GetXaxis()->SetTitle("Date & Time");
  gr1_time_start->GetXaxis()->SetTimeDisplay(kTRUE);
  gr1_time_start->GetXaxis()->SetTimeOffset(0,"gmt");    
  gr1_time_start->GetXaxis()->SetTimeFormat("%Y-%m-%d");  
  gr1_time_start->GetYaxis()->SetTitle("Parameter"); 
  gr1_time_start->Draw("AP");
  l1_time->Draw("same");

  c1_time->Write("45_210_nr_hr, time");
  
  TCanvas *c2_time = new TCanvas("c2_time","RDAnalysis_rp_45_210_fr_hr_time",700,500);
  c2_time->SetGrid();
  c2_time->SetLogy();
  TGraphErrors *gr2_time = new TGraphErrors(ncertainruns,&(start_times[0]),&(ycertainruns_fr_rp[0]),NULL,&(eycertainruns_fr_rp[0]));
  TF1 *l2_time = new TF1("ref","55.1898",start_times.front(),start_times.back());  
  gr2_time->SetMarkerStyle(21);
  gr2_time->SetMarkerSize(0.6);
  gr2_time->SetMaximum(65.);
  gr2_time->SetTitle("RDAnalysis_rp_45_210_fr_hr");
  gr2_time->GetXaxis()->SetTitle("Date & Time");
  gr2_time->GetXaxis()->SetTimeDisplay(kTRUE);
  gr2_time->GetXaxis()->SetTimeOffset(0,"gmt");  
  gr2_time->GetXaxis()->SetTimeFormat("%Y-%m-%d");  
  gr2_time->GetYaxis()->SetTitle("Parameter"); 
  gr2_time->Draw("AP");
  l2_time->Draw("same");

  c2_time->Write("45_210_fr_hr, time");
  
//   calculates parameter for the 56 nr rp with all the "good" runs
// cout << ncertainruns << endl;

  for (int i=0;i<ncertainruns;i++){
    stringstream ss;
    ss << "Runs/DQM_V0001_CTPPS_R000" << certainruns[i] << ".root"; //name of the root file
//     cout << ss.str() << endl;
    TFile *file = TFile::Open((ss.str()).c_str());
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/TrackingStrip/sector 56/station 210/nr_hr/track XY profile"; //path to the histogram
    for (int binx=11;binx<18;binx++){ //compute integral over the selected area
      for (int biny=49;biny<51;biny++){
	if (certainruns[i]<274244){
	ycertainruns_nr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
      	}  
	else{
	ycertainruns_nr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny);
	}
      }
    }
    
    double temp_sum = 0;
    
    for (int binx=23;binx<36;binx++){ //compute integral over the reference area
      for (int biny=32;biny<37;biny++){
	if (certainruns[i]<274244){
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
	}
	else{
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny);   
	}
      }
    }    
//     ycertainruns_nr_rp.at(i) /= (((TH2F*) file->Get(ss.str().c_str()))->GetEntries() - ycertainruns_nr_rp.at(i));
    eycertainruns_nr_rp.at(i) = ycertainruns_nr_rp.at(i) * sqrt(1./ycertainruns_nr_rp.at(i) + 1./temp_sum) / temp_sum;
    ycertainruns_nr_rp.at(i) /= temp_sum;
//     eycertainruns_nr_rp.at(i) = sqrt(ycertainruns_nr_rp.at(i))/luminosities[i];
//     ycertainruns_nr_rp.at(i) /= luminosities[i];
//     cout << ycertainruns_nr_rp.at(i) << endl;
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/events per BX";
    temp_sum_normalized_nr_rp.at(i) = temp_sum;
    interaction_BX_nr_rp.at(i) = ((TH2F*) file->Get(ss.str().c_str()))->GetEntries();
//     etemp_sum_normalized_nr_rp.at(i) = temp_sum_normalized_nr_rp.at(i)*sqrt(1./((TH2F*) file->Get(ss.str().c_str()))->GetEntries() + 1./temp_sum);
    ss.str(""); //resets string
  }
  
//   calculates parameter for the 56 fr rp with all the "good" runs
// cout << ncertainruns << endl;
  for (int i=0;i<ncertainruns;i++){
    stringstream ss;
    ss << "Runs/DQM_V0001_CTPPS_R000" << certainruns[i] << ".root"; //name of the root file
//     cout << ss.str() << endl;
    TFile *file = TFile::Open((ss.str()).c_str());
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/TrackingStrip/sector 56/station 210/fr_hr/track XY profile"; //path to the histogram
    for (int binx=12;binx<19;binx++){ //compute integral over the selected area
      for (int biny=41;biny<43;biny++){
	if (certainruns[i]<274244){
	ycertainruns_fr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);	  
	}
	else {
	ycertainruns_fr_rp.at(i) += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny);  	  
	}
      }
    }
    
    double temp_sum = 0;
    for (int binx=24;binx<37;binx++){ //compute integral over the reference area
      for (int biny=24;biny<29;biny++){
	if (certainruns[i]<274244){
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx-1,biny);
	}
	else{
	temp_sum += ((TH2F*) file->Get(ss.str().c_str()))->GetBinContent(binx,biny);
	}  
      }
    }
//     ycertainruns_fr_rp.at(i) /= (((TH2F*) file->Get(ss.str().c_str()))->GetEntries() - ycertainruns_fr_rp.at(i));
    eycertainruns_fr_rp.at(i) = ycertainruns_fr_rp.at(i) * sqrt(1./ycertainruns_fr_rp.at(i) + 1./temp_sum) / temp_sum;
    ycertainruns_fr_rp.at(i) /= temp_sum;
//     eycertainruns_fr_rp.at(i) = sqrt(ycertainruns_fr_rp.at(i))/luminosities[i];
//     ycertainruns_fr_rp.at(i) /= luminosities[i]; 
//     cout << ycertainruns_nr_rp.at(i) << endl;
    ss.str(""); //resets string
    ss << "DQMData/Run " << certainruns[i] << "/CTPPS/Run summary/events per BX";
    temp_sum_normalized_fr_rp.at(i) = temp_sum;
    interaction_BX_fr_rp.at(i) = ((TH2F*) file->Get(ss.str().c_str()))->GetEntries();
//     etemp_sum_normalized_fr_rp.at(i) = temp_sum_normalized_fr_rp.at(i)*sqrt(1./((TH2F*) file->Get(ss.str().c_str()))->GetEntries() + 1./temp_sum);
    ss.str(""); //resets string  
  }

	gDirectory = f_out;

// graphics part
//c3,c4 plots of parameter vs total integrated luminosity
//check_ref_area1 plots of correlation between hits in the reference area and # interactions per BX
//c3_time,c4_time plots of parameter vs time

  TCanvas *c3 = new TCanvas("c3","RDAnalysis_rp_56_210_nr_hr",700,500);
  c3->SetGrid();
  c3->SetLogy();
//   TGraphErrors *gr3 = new TGraphErrors(ncertainruns,certainruns,&(ycertainruns_nr_rp[0]),NULL,&(eycertainruns_nr_rp[0]));
  TGraphErrors *gr3 = new TGraphErrors(ncertainruns,&(luminosities_int[0]),&(ycertainruns_nr_rp[0]),NULL,&(eycertainruns_nr_rp[0]));
  TF1 *l3 = new TF1("ref","13.8614",0,7);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerSize(0.6);  
  gr3->SetTitle("RDAnalysis_rp_56_210_nr_hr");
  gr3->GetXaxis()->SetTitle("Integrated Luminosity (fb^{-1})");  
  gr3->GetYaxis()->SetTitle("Parameter"); 
  gr3->Draw("AP");  
  l3->Draw("same");

  c3->Write("56_210_nr_hr, lumi");
  
  TCanvas *c4 = new TCanvas("c4","RDAnalysis_rp_56_210_fr_hr",700,500);
  c4->SetGrid();
  c4->SetLogy();
//   TGraphErrors *gr4 = new TGraphErrors(ncertainruns,certainruns,&(ycertainruns_fr_rp[0]),NULL,&(eycertainruns_fr_rp[0]));   
  TGraphErrors *gr4 = new TGraphErrors(ncertainruns,&(luminosities_int[0]),&(ycertainruns_fr_rp[0]),NULL,&(eycertainruns_fr_rp[0]));   
  TF1 *l4 = new TF1("ref","15.1518",0,7);
  gr4->SetMarkerStyle(21);
  gr4->SetMarkerSize(0.6);
  gr4->SetTitle("RDAnalysis_rp_56_210_fr_hr");
  gr4->GetXaxis()->SetTitle("Integrated Luminosity (fb^{-1})");  
  gr4->GetYaxis()->SetTitle("Parameter"); 
  gr4->Draw("AP");
  l4->Draw("same");

  c4->Write("56_210_fr_hr, lumi");
  
//Graphics to check if the ratio between reference area hits and events per bunch crossing is constant
  TCanvas *check_ref_area1 = new TCanvas("check_ref_area1","rp_56_210");
  check_ref_area1->Divide(2);
  TPad* check_ref_area1_1 = (TPad*)(check_ref_area1->GetPrimitive("check_ref_area1_1"));
  TPad* check_ref_area1_2 = (TPad*)(check_ref_area1->GetPrimitive("check_ref_area1_2"));
  check_ref_area1_1->cd();
  TGraph *check_ref_area1_fr = new TGraph(ncertainruns,&(temp_sum_normalized_fr_rp[0]),&(interaction_BX_fr_rp[0]));   
  check_ref_area1_fr->SetMarkerStyle(21);
//   cout << check_ref_area1_fr->GetCorrelationFactor();
  check_ref_area1_fr->SetMarkerSize(0.6);
  check_ref_area1_fr->SetTitle("Reference area vs interactions correlation 56_210_fr_hr");
  check_ref_area1_fr->GetXaxis()->SetTitle("Hits in the reference area");
  check_ref_area1_fr->GetYaxis()->SetTitle("Total number of triggers");
  check_ref_area1_fr->GetYaxis()->SetTitleOffset(1.6);
  check_ref_area1_fr->Draw("AP");
  check_ref_area1_2->cd();
  TGraph *check_ref_area1_nr = new TGraph(ncertainruns,&(temp_sum_normalized_nr_rp[0]),&(interaction_BX_nr_rp[0]));   
  check_ref_area1_nr->SetMarkerStyle(21);
  check_ref_area1_nr->SetMarkerSize(0.6);
  check_ref_area1_nr->SetTitle("Reference area vs interactions correlation 56_210_nr_hr");
  check_ref_area1_nr->GetXaxis()->SetTitle("Hits in the reference area");
  check_ref_area1_nr->GetYaxis()->SetTitle("Total number of triggers");
  check_ref_area1_nr->GetYaxis()->SetTitleOffset(1.6);
  check_ref_area1_nr->Draw("AP");

  check_ref_area1->Write("56_210, check_ref_area");
  
  TCanvas *c3_time = new TCanvas("c3_time","RDAnalysis_rp_56_210_nr_hr_time",700,500);
  c3_time->SetGrid();
  c3_time->SetLogy();
  TGraphErrors *gr3_time = new TGraphErrors(ncertainruns,&(start_times[0]),&(ycertainruns_nr_rp[0]),NULL,&(eycertainruns_nr_rp[0]));
  TF1 *l3_time = new TF1("ref","13.8614",start_times.front(),start_times.back());
  gr3_time->SetMarkerStyle(21);
  gr3_time->SetMarkerSize(0.6);
  gr3_time->SetTitle("RDAnalysis_rp_56_210_nr_hr");
  gr3_time->GetXaxis()->SetTitle("Date & Time");
  gr3_time->GetXaxis()->SetTimeDisplay(kTRUE);
  gr3_time->GetXaxis()->SetTimeOffset(0,"gmt");  
  gr3_time->GetXaxis()->SetTimeFormat("%Y-%m-%d");  
  gr3_time->GetYaxis()->SetTitle("Parameter"); 
  gr3_time->Draw("AP");
  l3_time->Draw("same");

  c3_time->Write("56_210_nr_hr, time");
  
  TCanvas *c4_time = new TCanvas("c4_time","RDAnalysis_rp_56_210_fr_hr_time",700,500);
  c4_time->SetGrid();
  c4_time->SetLogy();
  TGraphErrors *gr4_time = new TGraphErrors(ncertainruns,&(start_times[0]),&(ycertainruns_fr_rp[0]),NULL,&(eycertainruns_fr_rp[0]));
  TF1 *l4_time = new TF1("ref","15.1518",start_times.front(),start_times.back());
  gr4_time->SetMarkerStyle(21);
  gr4_time->SetMarkerSize(0.6);
  gr4_time->SetTitle("RDAnalysis_rp_56_210_fr_hr");
  gr4_time->GetXaxis()->SetTitle("Date & Time");
  gr4_time->GetXaxis()->SetTimeDisplay(kTRUE);
  gr4_time->GetXaxis()->SetTimeOffset(0,"gmt");  
  gr4_time->GetXaxis()->SetTimeFormat("%Y-%m-%d");  
  gr4_time->GetYaxis()->SetTitle("Parameter"); 
  gr4_time->Draw("AP");
  l4_time->Draw("same");

  c4_time->Write("56_210_fr_hr, time");

  delete f_out;

  return 0;
}
