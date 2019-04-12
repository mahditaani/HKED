// C++ Includes
#include <iostream>
#include <algorithm>
#include <iomanip>

// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TFile.h"
#include "TApplication.h"
#include "TLatex.h"


// WCSim Includes
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

// Defines
#define PI 3.141592654
//#define MAXPMT 38448
//#define MAXPMTA 18604
#define WEIGHT 1
#define FILLW 0.0001
#define OFFSETID 25
#define OFFSETOD  25


#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>

//class TGRadioButton;

// Ensure that you have set the WCSIMDIR environment variable so llib can load the libraries required
/*
* To run this script just type:  root -l -x llib.C 'HKED.C("wcsim.root")'
* or replace wcsim.root your input filename
*/
std::string OutString(double n, int p = 2){


	// Create an output string stream
	std::ostringstream streamObj;

	// Set Fixed -Point Notation
	streamObj << std::fixed;

	// Set precision to 2 digits
	streamObj << std::setprecision(p);

	//Add double to stream
	streamObj << n;

	// Get string from output string stream
	std::string strObj = streamObj.str();
	return strObj;
}
std::string OutString(float n , int p = 2){
	return OutString((double)n,p);
}
std::string OutString(int n , int p = 0){
	return OutString((double)n, p);
}

bool checkPMT(int pmt, int low, int high) {

	if (pmt >= low && pmt <= high) return true;
	else return false;

}

class MyMainFrame {
   RQ_OBJECT("MyMainFrame")
private:
   TGMainFrame         *fMain;
   TRootEmbeddedCanvas *fEcanvasID;
	 TGHorizontalFrame *TopFrame;

	 TGTextButton *ButtonNext;
	 TGTextButton *ButtonPrev;
	 TGTextButton *ButtonSwitch;
	 TGTextButton *ButtonSave;
	 TGTextButton *ButtonExit;
	 TGTextEdit *TextBox;
	 std::string strings[14];
	 TFile *inFile;
	 int nEvent;
	 int ev =0;
	 TTree* geoTree;
	 WCSimRootGeom *geo;
	 TBranch *branchG;
	 TTree *wcsimTree;
	 TBranch *branch;
	 TBranch *branchOD;
	 int MAXPMT;
	 int MAXPMTA;
	 WCSimRootTrigger *wcsimTriggerID;
	 WCSimRootEvent *wcsimRootID;
	 WCSimRootTrigger *wcsimTriggerOD;
	 WCSimRootEvent *wcsimRootOD;
	 TH2Poly *blankID;
	 TH2Poly *blankOD;
	 TH2Poly *displayID;
	 TH2Poly *displayOD;
	 TCanvas *canvasID;
	 TCanvas *canvasOD;
	 std::string lineText;
	 double RadiusID;
	 double HeightID;
	 double RadiusOD;
	 double HeightOD;
	 bool idOn;
	 bool odOn;
	 bool bigID = true;
	 const int txtW = 20;
	 const int numW = 10;
	 int verbosity = 0;
	 const char *inFileName;
	 std::string fname;
	 double ene = 0;
	 double vtxx = 0;
	 double vtxy = 0;
	 double vtxz = 0;
	 double dirx = 0;
	 double diry = 0;
	 double dirz = 0;
	 int pmthitdigiid = 0;
	 int pmthitrawid = 0;
	 int hitdigiid = 0;
	 int hitrawid = 0;
	 int pmthitdigiod = 0;
	 int pmthitrawod = 0;
	 int hitdigiod = 0;
	 int hitrawod = 0;
	 int TopFrameHeight;
	 int TopFrameWidth;
	 int TopFrameX;
	 int TopFrameY;
	 int stringLength = 14;
public:
	MyMainFrame(std::string s);
   virtual ~MyMainFrame();
   void Prev();
	 void Next();
	 void Switch();
	 void Vision();
	 void Active();
	 void SaveCanvas();
	 void SetStrings();
	 void UpdateText();

};



void MyMainFrame::Active(){


	// Some nicely formatted text options
	std::cout << std::scientific; // This causes all numbers to be displayed in scientific notation.
	std::cout << std::setprecision(2); // Sets the decimal precision (no more than two decimal places)
	std::cout << std::left; // Sets the text justification to left
	const int txtW = 20; // Width of "box" holding text
	const int numW = 10; // Width of "box" holding numbers
	// Detector geometry details
	int MAXPMT = geo->GetWCNumPMT(); // Get the maximum number of PMTs in the ID
	int MAXPMTA = geo->GetODWCNumPMT(); // Get the maximum number of PMTs in the OD

	// Variables to read in from the root file ID
	float vtxX, vtxY, vtxZ; // vertex coordinate
	float dirX, dirY, dirZ; // particle momentum direction
	float energy; // particle energy
	int rawHitsID; // number of raw pmt hits
	int digiHitsID; // number of pmt digitised hits
	int numPMTsHitID; // Number of PMTs hit
	int numPMTsDigiHitID; // Number of PMTs with digihits
	// Variables to read in from the root file OD
	float vtxXOD, vtxYOD, vtxZOD; // vertex coordinate
	float dirXOD, dirYOD, dirZOD; // particle momentum direction
	float energyOD; // particle energy
	int rawHitsOD; // number of raw pmt hits
	int digiHitsOD; // number of pmt digitised hits
	int numPMTsHitOD; // Number of PMTs hit
	int numPMTsDigiHitOD; // Number of PMTs with digihits



	std::string command;
	wcsimTree->GetEntry(ev);
	wcsimTriggerID = wcsimRootID->GetTrigger(0);
	int numTriggersID = wcsimRootID->GetNumberOfEvents();
	int numSubTriggersID = wcsimRootID->GetNumberOfSubEvents();

	wcsimTriggerOD = wcsimRootOD->GetTrigger(0);
	int numTriggersOD = wcsimRootOD->GetNumberOfEvents();
	int numSubTriggersOD = wcsimRootOD->GetNumberOfSubEvents();

	//event = ev;

	if (verbosity) { // output the information of each event

		 std::cout << "======================================================" << std::endl;
		 std::cout << "******************************************************" << std::endl;
		 std::cout << "======================================================" << std::endl;
		 std::cout <<  std::left << std::setw(txtW) << "Event:" << std::right << std::setw(numW) << ev << std::endl;
		 std::cout <<  std::left << std::setw(txtW) << "Triggers:     " <<  std::right << std::setw(numW) << numTriggersID << std::endl;
		 std::cout <<  std::left << std::setw(txtW) << "Sub Triggers: " <<  std::right << std::setw(numW) << numSubTriggersID << std::endl;
		 std::cout << "======================================================" << std::endl;
		 std::cout <<  std::left << std::setw(txtW) << "EventOD:" << std::right << std::setw(numW) << ev << std::endl;
		 std::cout <<  std::left << std::setw(txtW) << "TriggersOD:     " <<  std::right << std::setw(numW) << numTriggersID << std::endl;
		 std::cout <<  std::left << std::setw(txtW) << "Sub TriggersOD: " <<  std::right << std::setw(numW) << numSubTriggersID << std::endl;

	} // End of if statement

	// Create copies of the blank event display
	displayID = (TH2Poly*) blankID->Clone();
	displayOD = (TH2Poly*) blankOD->Clone();

	for (int nTrig = 0; nTrig < numTriggersID; nTrig++){


		// ID
		wcsimTriggerID = wcsimRootID->GetTrigger(nTrig);
		int numTracksID = wcsimTriggerID->GetNtrack();
		WCSimRootTrack * trackID = (WCSimRootTrack*) wcsimTriggerID->GetTracks()->At(0);

		vtxX = wcsimTriggerID->GetVtx(0);
		vtxY = wcsimTriggerID->GetVtx(1);
		vtxZ = wcsimTriggerID->GetVtx(2);
		dirX = trackID->GetDir(0);
		dirY = trackID->GetDir(1);
		dirZ = trackID->GetDir(2);
		energy = trackID->GetE();
		rawHitsID = 0;
		digiHitsID = 0;

		numPMTsHitID = wcsimTriggerID->GetNcherenkovhits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
		numPMTsDigiHitID = wcsimTriggerID->GetNcherenkovdigihits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
		// END OF ID

		// ID
		wcsimTriggerOD = wcsimRootOD->GetTrigger(nTrig);
		int numTracksOD = wcsimTriggerOD->GetNtrack();
		WCSimRootTrack * trackOD = (WCSimRootTrack*) wcsimTriggerOD->GetTracks()->At(0);

		vtxXOD = wcsimTriggerOD->GetVtx(0);
		vtxYOD = wcsimTriggerOD->GetVtx(1);
		vtxZOD = wcsimTriggerOD->GetVtx(2);
		dirXOD = trackOD->GetDir(0);
		dirYOD = trackOD->GetDir(1);
		dirZOD = trackOD->GetDir(2);
		energyOD = trackOD->GetE();
		rawHitsOD = 0;
		digiHitsOD = 0;

		numPMTsHitOD = wcsimTriggerOD->GetNcherenkovhits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
		numPMTsDigiHitOD = wcsimTriggerOD->GetNcherenkovdigihits(); //Returns the number of PMTs with a true hit (photon or dark noise) (QE applied)
		// END OF ID



		// Work out the number of photons that hit the PMTs
		for (int i = 0; i < numPMTsHitID; i++){

		WCSimRootCherenkovHit *cherenkovHitID = (WCSimRootCherenkovHit*) wcsimTriggerID->GetCherenkovHits()->At(i);
		float tmpRawHitsID = cherenkovHitID->GetTotalPe(1);
		rawHitsID += tmpRawHitsID;


		} // End of for loop working out number of photons in PMTs
		// Work out the number of photons that hit the PMTs
		for (int i = 0; i < numPMTsHitOD; i++){

		WCSimRootCherenkovHit *cherenkovHitOD = (WCSimRootCherenkovHit*) wcsimTriggerOD->GetCherenkovHits()->At(i);
		float tmpRawHitsOD = cherenkovHitOD->GetTotalPe(1);
		rawHitsOD += tmpRawHitsOD;


		} // End of for loop working out number of photons in PMTs

		// Work out the number of digitised hits in the PMTs
		for (int i = 0; i < numPMTsDigiHitID; i++){

		WCSimRootCherenkovDigiHit *cherenkovDigiHitID = (WCSimRootCherenkovDigiHit*) wcsimTriggerID->GetCherenkovDigiHits()->At(i);
		float tmpDigiHitsID = cherenkovDigiHitID->GetQ();

		digiHitsID += tmpDigiHitsID;
		// Find the tube and work out its position
		int tubeID = cherenkovDigiHitID->GetTubeId() -1;
		double tube[3];
		int cylLoc = geo->GetPMT(tubeID).GetCylLoc();
		tube[0] = geo->GetPMT(tubeID).GetPosition(0);
		tube[1] = geo->GetPMT(tubeID).GetPosition(1);
		tube[2] = geo->GetPMT(tubeID).GetPosition(2);

		//Top ID
		if ( cylLoc == 0){
			displayID->Fill(tube[0],tube[1] + RadiusID + HeightID/2,tmpDigiHitsID);
		}
		//Bot ID
		else if ( cylLoc == 2){
			displayID->Fill(tube[0],-(HeightID/2 +RadiusID +tube[1]),tmpDigiHitsID);
		}
		//Barrel ID
		else {
			double l = sqrt( pow((0 - tube[0]),2) + pow((-RadiusID - tube[1]),2));
			double angle = 2*asin(l/(2*RadiusID));
			double length = angle*RadiusID ;
			if (tube[0]<0) length *= -1;
			displayID->Fill( length, tube[2],tmpDigiHitsID);
		}


		} // End of for loop working out number of digitized hits in PMTs

		// Work out the number of digitised hits in the PMTs
		for (int i = 0; i < numPMTsDigiHitOD; i++){

		WCSimRootCherenkovDigiHit *cherenkovDigiHitOD = (WCSimRootCherenkovDigiHit*) wcsimTriggerOD->GetCherenkovDigiHits()->At(i);
		float tmpDigiHitsOD = cherenkovDigiHitOD->GetQ();

		digiHitsOD += tmpDigiHitsOD;
		// Find the tube and work out its position
		int tubeOD = MAXPMT + cherenkovDigiHitOD->GetTubeId() -1; /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double tube[3];
		int cylLoc = geo->GetPMT(tubeOD).GetCylLoc();
		tube[0] = geo->GetPMT(tubeOD).GetPosition(0);
		tube[1] = geo->GetPMT(tubeOD).GetPosition(1);
		tube[2] = geo->GetPMT(tubeOD).GetPosition(2);

		//Top OD
		if ( cylLoc == 5){
			displayOD->Fill(tube[0],tube[1] + RadiusOD + HeightOD/2 + 100,tmpDigiHitsOD);
		}
		//Bot OD
		else if ( cylLoc == 3){
			displayOD->Fill(tube[0],-(HeightOD/2 +RadiusOD +tube[1] + 100 ),tmpDigiHitsOD);
		}
		//Barrel OD
		else {

			double l = sqrt( pow((0 - tube[0]),2) + pow((-RadiusOD - tube[1]),2));
			double angle = 2*asin(l/(2*RadiusOD));
			double length = angle*RadiusOD ;
			if (tube[0]<0) length *= -1;
			displayOD->Fill( length, tube[2],tmpDigiHitsOD);
		}



		} // End of for loop working out number of digitized hits in PMTs
		if (verbosity) { // output the information of each event
			std::cout << "======================================================" << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Event:" <<  std::right << std::setw(numW) << ev << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "VTX (x y z)[cm]:" <<  std::right << std::setw(numW) << vtxX << " " <<  std::setw(numW) << vtxY << " " <<  std::setw(numW) << vtxZ << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Dir (x y z)[norm]:" <<  std::right << std::setw(numW) << dirX << " " <<  std::setw(numW) << dirY << " " <<  std::setw(numW) << dirZ << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Energy [MeV]:" <<  std::right << std::setw(numW) << energy << std::endl;
			std::cout << "--------------- ID ---------------------" << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "PMTs Hit Digi:" <<  std::right << std::setw(numW) << numPMTsDigiHitID << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "PMTs Hit Raw:" <<  std::right << std::setw(numW) << numPMTsHitID << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Raw Hits:" <<  std::right << std::setw(numW) << rawHitsID << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Digi Hits:" <<  std::right << std::setw(numW) << digiHitsID << std::endl;
			std::cout << "--------------- OD ---------------------" << std::endl;
			std::cout << std::left ; // change text justification back to left.
			std::cout <<  std::left << std::setw(txtW) << "PMTs Hit Digi:" <<  std::right << std::setw(numW) << numPMTsDigiHitOD << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "PMTs Hit Raw:" <<  std::right << std::setw(numW) << numPMTsHitOD << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Raw Hits:" <<  std::right << std::setw(numW) << rawHitsOD << std::endl;
			std::cout <<  std::left << std::setw(txtW) << "Digi Hits:" <<  std::right << std::setw(numW) << digiHitsOD << std::endl;
			std::cout << std::left ; // change text justification back to left.
		} // End of if statement


	} // End of loop over triggers


	if (idOn || odOn){

		canvasID->cd();

		TPad *bigPad = new TPad("big", "", 0,0,1, 1);
		bigPad->SetFillColor(kBlack);
		bigPad->SetRightMargin(0.00);
		bigPad->SetLeftMargin(0.);
		bigPad->SetTopMargin(0.00);
		bigPad->SetBottomMargin(0.);

		bigPad->Draw();
		bigPad->cd();

		if( !bigID && odOn){
			displayOD->Draw("COLZ");
		} else if( bigID && idOn){
			displayID->Draw("COLZ");
		}

		for(int gh = 0; gh < 4; gh++){

			float xtxt = 0.025;
			float ytxt = 0.8 - gh*0.03;
			std::string tmpString = "#color[2]{" + strings[gh] + "}";
			TLatex *textInfo = new TLatex(xtxt,ytxt, tmpString.c_str());
			textInfo->SetLineWidth(2);
			textInfo->SetNDC(kTRUE);
			textInfo->SetTextSize(0.02);
			textInfo->Draw();

		}

		for(int gh = 4; gh < 14; gh++){

			float xtxt = 0.025;
			float ytxt = 0.34 - (gh-4)*0.03;
			std::string tmpString = "#color[5]{" + strings[gh] + "}";
			TLatex *textInfo = new TLatex(xtxt,ytxt, tmpString.c_str());
			textInfo->SetLineWidth(2);
			textInfo->SetNDC(kTRUE);
			textInfo->SetTextSize(0.02);
			textInfo->Draw();

		}


		bigPad->Modified();
		canvasID->cd();

		TPad *smallPad = new TPad("small", "", 0.68,0.65,1, 1);
		smallPad->SetFillColor(kBlack);
		smallPad->SetRightMargin(0.00);
		smallPad->SetLeftMargin(0.);
		smallPad->SetTopMargin(0.00);
		smallPad->SetBottomMargin(0.);

		smallPad->Draw();
		smallPad->cd();
		if( bigID && odOn){
			displayOD->Draw("COLZ");
		} else if( !bigID && idOn){
			displayID->Draw("COLZ");
		}

		smallPad->Modified();
		canvasID->cd();

		canvasID->Modified();
		canvasID->cd();
		canvasID->SetSelected(canvasID);
		canvasID->Update();

	}


	ene = energy;
	vtxx = vtxX;
	vtxy = vtxY;
	vtxz = vtxZ;
	dirx = dirX;
	diry = dirY;
	dirz = dirZ;
	pmthitdigiid = numPMTsDigiHitID;
	pmthitrawid = numPMTsHitID;
	hitdigiid = digiHitsID;
	hitrawid = rawHitsID;
	pmthitdigiod = numPMTsDigiHitOD;
	pmthitrawod = numPMTsHitOD;
	hitdigiod = digiHitsOD;
	hitrawod = rawHitsOD;

}

void MyMainFrame::Vision(){

	gStyle->SetOptStat(0); // Remove stats from our histograms
	gStyle->SetPalette(1); // Use the rainbow colour palette
	// Some nicely formatted text options
	std::cout << std::scientific; // This causes all numbers to be displayed in scientific notation.
	std::cout << std::setprecision(2); // Sets the decimal precision (no more than two decimal places)
	std::cout << std::left; // Sets the text justification to left


	ev = 0;

	// Open the WCSim file
	inFile = new TFile(inFileName, "READ");
	if ( !inFile->IsOpen() ){
		std::cout << "Error: could not open input file \"" << inFileName << "\"." <<std::endl;

	} else if (verbosity) {
		std::cout << "Input file: " << inFileName << std::endl;
	}



	// Get a pointer to the tree from the input file
	wcsimTree = (TTree*) inFile->Get("wcsimT");

	// Get the number of events in the tree
	nEvent = wcsimTree->GetEntries();
	if (verbosity) { std::cout << "Number of events: "<< nEvent << std::endl;}

	// Create a WCSimRootEvent to put stuff from the tree in
	wcsimRootID = new WCSimRootEvent();
	wcsimRootOD = new WCSimRootEvent();

	// Set the branch address for reading from the tree
	branch = wcsimTree->GetBranch("wcsimrootevent");
	branchOD = wcsimTree->GetBranch("wcsimrootevent_OD");
	branch->SetAddress(&wcsimRootID);
	branchOD->SetAddress(&wcsimRootOD);

	// Force deletion to prevent memory leak
	wcsimTree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
	wcsimTree->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);

	// Load the geometry tree (only 1 "event")
	geoTree = (TTree*) inFile->Get("wcsimGeoT");
	geo = new WCSimRootGeom();
	branchG = geoTree->GetBranch("wcsimrootgeom");
	branchG->SetAddress(&geo);

	if (verbosity) {std::cout << "Geotree has " << geoTree->GetEntries() << " entries." << std::endl;}

	geoTree->GetEntry(0);


	canvasID = fEcanvasID->GetCanvas();

	canvasID->SetRightMargin(0.00);
	canvasID->SetLeftMargin(0.);
	canvasID->SetTopMargin(0.00);
	canvasID->SetBottomMargin(0.);

	// Plot parameters
	int nBinID = 100; // Granularity for ID plot
	int nBinOD = 100; // Granularity for OD plot
	double dimX = 10000;  // 2D histogram half length of x axis
	double dimY = 10000;  // 2D histogram half length of y axis
	double dimZ = 10000; // 2D histogram half length of z axis

	// Detector geometry details
	MAXPMT = geo->GetWCNumPMT(); // Get the maximum number of PMTs in the ID
	MAXPMTA = geo->GetODWCNumPMT(); // Get the maximum number of PMTs in the OD

	idOn = true; // Boolean to keep track of whether or not the ID was constructed (sometimes turned off for speed) true = on , false = off
	odOn = true; // Boolean to keep track of whether or not the ID was constructed (sometimes turned off for speed) true = on , false = off
	if (MAXPMT == 0 ) {idOn = false;}
	if (MAXPMTA == 0 ) {odOn = false;}

	RadiusID = 0;
	RadiusOD = 0;
	HeightID = 0;
	HeightOD = 0;


	displayID = new TH2Poly ("displayID", "displayID", -dimX, dimX, -dimZ, dimZ);
	displayOD = new TH2Poly ("displayOD", "displayOD", -dimX, dimX, -dimZ, dimZ);
	blankID = new TH2Poly ("blankID", "", -dimX, dimX, -dimZ, dimZ);
	blankOD = new TH2Poly ("blankOD", "", -dimX, dimX, -dimZ, dimZ);

	// Remove Axis Lables
	blankID->GetXaxis()->SetLabelOffset(999);
  blankID->GetXaxis()->SetLabelSize(0);
	blankOD->GetXaxis()->SetLabelOffset(999);
	blankOD->GetXaxis()->SetLabelSize(0);
	blankID->GetYaxis()->SetLabelOffset(999);
	blankID->GetYaxis()->SetLabelSize(0);
	blankOD->GetYaxis()->SetLabelOffset(999);
	blankOD->GetYaxis()->SetLabelSize(0);



	// Find a barrel pmt
	bool barrelID = false; // Boolean to see if a barrel ID pmt has been found
	bool barrelOD = false;  // Boolean to see if a barrel OD pmt has been found
	bool capID = false;  // Boolean to see if a cap ID pmt has been found
	bool capOD = false;  // Boolean to see if a cap OD pmt has been found
	int barrelPMTID = -1; // PMT number of the barrel ID PMT
	int barrelPMTOD = -1; // PMT number of the barrel OD PMT
	int capPMTID = -1; // PMT number of the cap ID PMT
	int capPMTOD = -1; // PMT number of the cap OD PMT
	int pmtCount = 0; // Number used to count through all of the PMTs

	if (!idOn) {barrelID = true; capID = true;} // If no ID constructed, don't look for ID PMTs
	if (!odOn) {barrelOD = true; capOD = true;} // If no OD constructed, don't look for OD PMTs

	while (!barrelID || !barrelOD || !capID || !capOD ){ // Loop to look for barrel and cap PMTs to work out the radius and height respectively.

		//std::cout <<"Looking for barrel ID PMT: " << geo->GetPMT(pmtCount).GetCylLoc()<<std::endl;
		if ( !barrelID && ( geo->GetPMT(pmtCount).GetCylLoc() == 1)  ) {barrelID = true; barrelPMTID = pmtCount; }
		if ( !barrelOD && (geo->GetPMT(pmtCount).GetCylLoc() == 4 )  ) {barrelOD = true; barrelPMTOD = pmtCount; }
		if ( !capID && (geo->GetPMT(pmtCount).GetCylLoc() == 0 || geo->GetPMT(pmtCount).GetCylLoc() == 2)  ) {capID = true; capPMTID = pmtCount; }
		if ( !capOD && (geo->GetPMT(pmtCount).GetCylLoc() == 3 || geo->GetPMT(pmtCount).GetCylLoc() == 5)  ) {capOD = true; capPMTOD = pmtCount; }
		pmtCount++;
		//pmtCount+= 10; // Can speed up this process by checking PMTs in multiples higher than 1
	}


	if (idOn) { // If ID is on, check the PMTs are correct and set the height and radius.
		if (checkPMT(barrelPMTID, 0, MAXPMT -1 ) && checkPMT(capPMTID, 0, MAXPMT -1) ) {
			// Set the radius and height of the ID using the PMTs' positions.
			RadiusID = sqrt( pow( geo->GetPMT(barrelPMTID).GetPosition(0),2) + pow(geo->GetPMT(barrelPMTID).GetPosition(1),2) );
			HeightID = 2*(abs(geo->GetPMT(capPMTID).GetPosition(2)));

		}
		else {
			std::cerr << "Can not understand the tank geometry. Exiting..." << std::endl;
			exit(1);
		}
	}

	if (odOn) { // If OD is on, check the PMTs are correct and set the height and radius.
		if (checkPMT(barrelPMTOD, MAXPMT, MAXPMT + MAXPMTA-1) && checkPMT(capPMTOD, MAXPMT, MAXPMT + MAXPMTA - 1) ) {
			// Set the radius and height of the ID and OD using the PMTs' positions.
			RadiusOD = sqrt( pow( geo->GetPMT(barrelPMTOD).GetPosition(0),2) + pow(geo->GetPMT(barrelPMTOD).GetPosition(1),2) );
			HeightOD = 2*(abs(geo->GetPMT(capPMTOD).GetPosition(2)));

		}
		else {
			std::cerr << "Can not understand the tank geometry. Exiting..." << std::endl;
			exit(1);
		}
	}


	if (idOn) {
		std::cout << "Barrel Radius (ID) is: " << RadiusID <<std::endl;
		std::cout << "Barrel Height (ID) " <<  HeightID <<std::endl;
	}
	if (odOn) {
		std::cout << "Barrel Radius (OD) is: " << RadiusOD <<std::endl;
		std::cout << "Barrel Height (OD) " <<  HeightOD <<std::endl;
	}

	// Fill up the cylinder shape for OD hits
	for ( int i = MAXPMT; i < MAXPMT + MAXPMTA; i++){
		double tbe[3];
		int cylLoc = geo->GetPMT(i).GetCylLoc();
		tbe[0] = geo->GetPMT(i).GetPosition(0);
		tbe[1] = geo->GetPMT(i).GetPosition(1);
		tbe[2] = geo->GetPMT(i).GetPosition(2);

		//Top
		if ( cylLoc == 5){
			blankOD->AddBin(tbe[0] - OFFSETOD, tbe[1] + RadiusOD + HeightOD/2 - OFFSETOD,
			tbe[0] + OFFSETOD, tbe[1] + RadiusOD + HeightOD/2 + OFFSETOD);
		}
		//Bot
		else if ( cylLoc == 3){
			blankOD->AddBin(tbe[0] - OFFSETOD, -HeightOD/2 - RadiusOD -tbe[1] - OFFSETOD,
				tbe[0] + OFFSETOD, -HeightOD/2 - RadiusOD -tbe[1] + OFFSETOD);
		}
		//Barrel
		else {
			double l = sqrt( pow((0 - tbe[0]),2) + pow((-RadiusOD - tbe[1]),2));
			double angle = 2*asin(l/(2*RadiusOD));
			double length = angle*RadiusOD ;
			if (tbe[0]<0) length *= -1;
			blankOD->AddBin( length - OFFSETOD, tbe[2] - OFFSETOD, length + OFFSETID, tbe[2] + OFFSETID );
		}

	} // End of for loop filling OD PMT hits

	// Fill up the cylinder shape for ID hits
	for ( int i = 0; i < MAXPMT; i++){
		double tbe[3];
		int cylLoc = geo->GetPMT(i).GetCylLoc();
		tbe[0] = geo->GetPMT(i).GetPosition(0);
		tbe[1] = geo->GetPMT(i).GetPosition(1);
		tbe[2] = geo->GetPMT(i).GetPosition(2);

		//Top
		if ( cylLoc == 0){
			blankID->AddBin(tbe[0] - OFFSETID, tbe[1] + RadiusID + HeightID/2 - OFFSETID,
				tbe[0] + OFFSETID, tbe[1] + RadiusID + HeightID/2 + OFFSETID);
		}
		//Bot
		else if ( cylLoc == 2){
			blankID->AddBin(tbe[0] - OFFSETID, -(HeightID/2 + RadiusID + tbe[1]) - OFFSETID,
			tbe[0] + OFFSETID, -(HeightID/2 + RadiusID + tbe[1]) + OFFSETID);
		}
		//Barrel
		else {
			double l = sqrt( pow((0 - tbe[0]),2) + pow((-RadiusID - tbe[1]),2));
			double angle = 2*asin(l/(2*RadiusID));
			double length = angle*RadiusID ;
			if (tbe[0]<0) length *= -1;
			blankID->AddBin( length- OFFSETID, tbe[2] - OFFSETID, length + OFFSETID, tbe[2] + OFFSETID);
		}


	} // End of for loop filling ID PMT hits


}

MyMainFrame::MyMainFrame(string s) {
	 inFileName = s.c_str();

   // Create a main frame
	 fMain = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
	 fMain->SetName("EventDisplay");
	 fMain->SetLayoutBroken(kTRUE);
	 fMain->SetEditable(kFALSE);

	 // top frame
   TopFrame = new TGHorizontalFrame(fMain,608,432,kHorizontalFrame);
   TopFrame->SetName("TopFrame");
   TopFrame->SetLayoutBroken(kTRUE);

   fMain->AddFrame(TopFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   TopFrame->MoveResize(0,0,1500,900);

	 TopFrameX = TopFrame->GetX();
	 TopFrameY = TopFrame->GetY();
	 TopFrameWidth = TopFrame->GetWidth();
	 TopFrameHeight = TopFrame->GetHeight();

   // Create canvas widget to go in the top frame
   fEcanvasID = new TRootEmbeddedCanvas("EcanvasID",TopFrame,100,100);
	 fEcanvasID->SetAutoFit();
	 fEcanvasID->MoveResize(TopFrameX, TopFrameY, TopFrameWidth, TopFrameHeight);
   TopFrame->AddFrame(fEcanvasID, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));

	 fEcanvasID->MoveResize(TopFrameX, TopFrameY, TopFrameWidth, TopFrameHeight);

   ButtonNext = new TGTextButton(fMain, "Next",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
	 ButtonNext->SetToolTipText("Displays the next event");
   ButtonNext->SetTextJustify(36);
   ButtonNext->SetMargins(0,0,0,0);
   ButtonNext->SetWrapLength(-1);
   fMain->AddFrame(ButtonNext, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   ButtonNext->MoveResize(144,920,92,24);

   ButtonPrev = new TGTextButton(fMain,"Previous",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
	 ButtonPrev->SetToolTipText("Displays the previous event");
   ButtonPrev->SetTextJustify(36);
   ButtonPrev->SetMargins(0,0,0,0);
   ButtonPrev->SetWrapLength(-1);
   fMain->AddFrame(ButtonPrev, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   ButtonPrev->MoveResize(32,920,92,24);

	 ButtonSwitch = new TGTextButton(fMain,"Switch",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
	 ButtonSwitch->SetToolTipText("Switches the order of the ID and OD.");
   ButtonSwitch->SetTextJustify(36);
   ButtonSwitch->SetMargins(0,0,0,0);
   ButtonSwitch->SetWrapLength(-1);
   fMain->AddFrame(ButtonSwitch, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   ButtonSwitch->MoveResize(256,920,92,24);

   ButtonSave = new TGTextButton(fMain,"Save",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
	 ButtonSave->SetToolTipText("Saves the displays, as you see them, to png files");
   ButtonSave->SetTextJustify(36);
   ButtonSave->SetMargins(0,0,0,0);
   ButtonSave->SetWrapLength(-1);
   fMain->AddFrame(ButtonSave, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   ButtonSave->MoveResize(368,920,92,24);

	 ButtonExit = new TGTextButton(fMain,"Exit","gApplication->Terminate(0)",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   ButtonExit->SetTextJustify(36);
   ButtonExit->SetMargins(0,0,0,0);
   ButtonExit->SetWrapLength(-1);

   fMain->AddFrame(ButtonExit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   ButtonExit->MoveResize(480,920,92,24);

   fMain->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fMain->MapSubwindows();

   fMain->Resize(fMain->GetDefaultSize());
   fMain->MapWindow();
   fMain->Resize(1500,950);



	 // Button Commands
   ButtonPrev->Connect("Clicked()","MyMainFrame",this,"Prev()");
   ButtonNext->Connect("Clicked()","MyMainFrame",this,"Next()");
	 ButtonSwitch->Connect("Clicked()","MyMainFrame",this,"Switch()");
	 ButtonSave->Connect("Clicked()","MyMainFrame",this,"SaveCanvas()");

	 SetStrings(); // Format the text to be seen in the bottom left

   // Set a name to the main frame
   fMain->SetWindowName("HK Event Display");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());


   // Map main frame
   fMain->MapWindow();

	 // Initialize other variables and plots
	 Vision();
	 Active();
	 UpdateText();
	 Active();

}

void MyMainFrame::SetStrings(){


	strings[0] = "Event";
	strings[0].insert(strings[0].end(), stringLength - strings[0].length(), '.');
	strings[0] = strings[0] + ": " + OutString(ev+1) + "/" + OutString(nEvent);

	strings[1] = "Energy";
	strings[1].insert(strings[1].end(), stringLength - strings[1].length(), '.');
	strings[1] = strings[1] + ": " + OutString(ene) + " MeV";

	strings[2] = "vtx [cm]";
	strings[2].insert(strings[2].end(), stringLength - strings[2].length(), '.');
	strings[2] = strings[2] + ": " + OutString(vtxx) + ", " + OutString(vtxy) + ", " + OutString(vtxz);

	strings[3] = "dir [norm]";
	strings[3].insert(strings[3].end(), stringLength - strings[3].length(), '.');
	strings[3] = strings[3] + ": " + OutString(dirx) + ", " + OutString(diry) + ", " + OutString(dirz);

	strings[4] = "---------ID-------------";

	strings[5] = "PMTs Hit Digi\t";
	strings[5].insert(strings[5].end(), stringLength - strings[5].length(), '.');
	strings[5] = strings[5] + ": " + OutString(pmthitdigiid);

	strings[6] = "PMTs Hit Raw";
	strings[6].insert(strings[6].end(), stringLength - strings[6].length(), '.');
	strings[6] = strings[6] + ": " + OutString(pmthitrawid);

	strings[7] = "Digi Hits";
	strings[7].insert(strings[7].end(), stringLength - strings[7].length(), '.');
	strings[7] = strings[7] + ": " + OutString(hitdigiid);

	strings[8] = "Raw Hits";
	strings[8].insert(strings[8].end(), stringLength - strings[8].length(), '.');
	strings[8] = strings[8] + ": " + OutString(hitrawid);

	strings[9] = "---------OD-------------";

	strings[10] = "PMTs Hit Digi";
	strings[10].insert(strings[10].end(), stringLength - strings[10].length(), '.');
	strings[10] = strings[10] + ": " + OutString(pmthitdigiod);

	strings[11] = "PMTs Hit Raw";
	strings[11].insert(strings[11].end(), stringLength - strings[11].length(), '.');
	strings[11] = strings[11] + ": " + OutString(pmthitrawod);

	strings[12] = "Digi Hits";
	strings[12].insert(strings[12].end(), stringLength - strings[12].length(), '.');
	strings[12] = strings[12] + ": " + OutString(hitdigiod);

	strings[13] = "Raw Hits";
	strings[13].insert(strings[13].end(), stringLength - strings[13].length(), '.');
	strings[13] = strings[13] + ": " + OutString(hitrawod);


}
void MyMainFrame::UpdateText(){

	wcsimTriggerID = wcsimRootID->GetTrigger(0);
	WCSimRootTrack * trackID = (WCSimRootTrack*) wcsimTriggerID->GetTracks()->At(0);
	vtxx = wcsimTriggerID->GetVtx(0);
	vtxy = wcsimTriggerID->GetVtx(1);
	vtxz = wcsimTriggerID->GetVtx(2);
 	dirx = trackID->GetDir(0);
 	diry = trackID->GetDir(1);
 	dirz = trackID->GetDir(2);
 	ene = trackID->GetE();
 	pmthitrawid = wcsimTriggerID->GetNcherenkovhits();
 	pmthitdigiid = wcsimTriggerID->GetNcherenkovdigihits();

 	wcsimTriggerOD = wcsimRootOD->GetTrigger(0);
 	pmthitrawod = wcsimTriggerOD->GetNcherenkovhits();
 	pmthitdigiod = wcsimTriggerOD->GetNcherenkovdigihits();
	int hitnum = 0;

	hitnum = wcsimTriggerOD->GetCherenkovDigiHits()->GetEntries();
	hitdigiod = 0;
	for (int jh = 0; jh < hitnum; jh++){
		WCSimRootCherenkovDigiHit *cherenkovDigiHitOD = (WCSimRootCherenkovDigiHit*) wcsimTriggerOD->GetCherenkovDigiHits()->At(jh);
		hitdigiod += cherenkovDigiHitOD->GetQ();
	}

	hitnum = wcsimTriggerOD->GetCherenkovHits()->GetEntries();
	hitrawod = 0;
	for (int jh = 0; jh < hitnum; jh++){
		WCSimRootCherenkovHit *cherenkovHitOD = (WCSimRootCherenkovHit*) wcsimTriggerOD->GetCherenkovHits()->At(jh);
		hitrawod += cherenkovHitOD->GetTotalPe(1);
	}

	hitnum = wcsimTriggerID->GetCherenkovDigiHits()->GetEntries();
	hitdigiid = 0;
	for (int jh = 0; jh < hitnum; jh++){
		WCSimRootCherenkovDigiHit *cherenkovDigiHitID = (WCSimRootCherenkovDigiHit*) wcsimTriggerID->GetCherenkovDigiHits()->At(jh);
		hitdigiid += cherenkovDigiHitID->GetQ();
	}

	hitnum = wcsimTriggerID->GetCherenkovHits()->GetEntries();
	hitrawid = 0;
	for (int jh = 0; jh < hitnum; jh++){
		WCSimRootCherenkovHit *cherenkovHitID = (WCSimRootCherenkovHit*) wcsimTriggerID->GetCherenkovHits()->At(jh);
		hitrawid += cherenkovHitID->GetTotalPe(1);
	}

 	SetStrings();
}


void MyMainFrame::Prev() {
	if (ev > 0){
		ev--;
		UpdateText();
		Active();
	}
}
void MyMainFrame::Next() {
	if (ev < nEvent - 1 ){
		ev++;
		UpdateText();
		Active();


	}
}
void MyMainFrame::Switch() {
		bigID = !bigID;
		UpdateText();
		Active();

}
void MyMainFrame::SaveCanvas(){
	if (idOn) {
		fEcanvasID->GetCanvas()->SaveAs("id.png");
	}
}

MyMainFrame::~MyMainFrame() {
   // Clean up used widgets: frames, buttons, layout hints
   fMain->Cleanup();
   delete fMain;
}
void HKED(std::string file) {

   // Popup the GUI...
   //new MyMainFrame(file, gClient->GetRoot(), 800, 600);
	 new MyMainFrame(file);

}
