// C++ Includes
#include <iostream>
#include <algorithm>
#include <iomanip>

// ROOT Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TApplication.h"


// WCSim Includes
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

// Defines
#define PI 3.141592654
//#define MAXPMT 38448
//#define MAXPMTA 18604
#define WEIGHT 1
#define FILLW 0.0001


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
	return OutString((double)n);
}
std::string OutString(int n , int p = 2){
	return OutString((double)n);
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
	 TRootEmbeddedCanvas *fEcanvasOD;
	 TGHorizontalFrame *frameLeft;
	 TGHorizontalFrame *frameRight;
	 TGVerticalFrame *EventDetailsFrame;
	 TGVerticalFrame *ButtonsFrame;
	 TGLabel *text[14];
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
	 TH2D *blankID;
	 TH2D *blankOD;
	 TH2D *displayID;
	 TH2D *displayOD;
	 TCanvas *canvasID;
	 TCanvas *canvasOD;
	 double RadiusID;
	 double HeightID;
	 double RadiusOD;
	 double HeightOD;
	 bool idOn;
	 bool odOn;
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
public:
   MyMainFrame(std::string s, const TGWindow *p,UInt_t w,UInt_t h);
   virtual ~MyMainFrame();
   void Prev();
	 void Next();
	 void Vision();
	 void Active();
	 void IDOnly();
	 void ODOnly();
	 void IDOD();
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
	displayID = (TH2D*) blankID->Clone();
	displayOD = (TH2D*) blankOD->Clone();

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

	if (odOn){ canvasOD->cd(); displayOD->Draw("COLZ"); canvasOD->Update(); }
	if (idOn){ canvasID->cd(); displayID->Draw("COLZ"); canvasID->Update(); }


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
	canvasOD = fEcanvasOD->GetCanvas();


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


	displayID = new TH2D ("displayID", "displayID", nBinID, -dimX, dimX, nBinID, -dimZ, dimZ);
	displayOD = new TH2D ("displayOD", "displayOD", nBinOD, -dimX, dimX, nBinOD, -dimZ, dimZ);
	blankID = new TH2D ("blankID", "ID", nBinID, -dimX, dimX, nBinID, -dimZ, dimZ);
	blankOD = new TH2D ("blankOD", "OD", nBinOD, -dimX, dimX, nBinOD, -dimZ, dimZ);


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
			blankOD->Fill(tbe[0],tbe[1] + RadiusOD + HeightOD/2 ,FILLW);
		}
		//Bot
		else if ( cylLoc == 3){
			blankOD->Fill(tbe[0],-HeightOD/2 - RadiusOD -tbe[1],FILLW);
		}
		//Barrel
		else {
			double l = sqrt( pow((0 - tbe[0]),2) + pow((-RadiusOD - tbe[1]),2));
			double angle = 2*asin(l/(2*RadiusOD));
			double length = angle*RadiusOD ;
			if (tbe[0]<0) length *= -1;
			blankOD->Fill( length, tbe[2],FILLW);
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
			blankID->Fill(tbe[0],tbe[1] + RadiusID + HeightID/2 ,FILLW);
		}
		//Bot
		else if ( cylLoc == 2){
			blankID->Fill(tbe[0],-(HeightID/2 + RadiusID + tbe[1]),FILLW);
		}
		//Barrel
		else {
			double l = sqrt( pow((0 - tbe[0]),2) + pow((-RadiusID - tbe[1]),2));
			double angle = 2*asin(l/(2*RadiusID));
			double length = angle*RadiusID ;
			if (tbe[0]<0) length *= -1;
			blankID->Fill( length, tbe[2],FILLW);
		}


	} // End of for loop filling ID PMT hits


}




MyMainFrame::MyMainFrame(string s,const TGWindow *p,UInt_t w,UInt_t h) {

	 inFileName = s.c_str();
   // Create a main frame
   fMain = new TGMainFrame(p,w,h);
	 TQObject::Connect("MyMainFrame", "CloseWindow()", 0, 0, "gApplication->Terminate(0)");

	 frameLeft = new TGHorizontalFrame(fMain,200,200);
	 frameRight = new TGHorizontalFrame(fMain,200,200);
	 EventDetailsFrame = new TGVerticalFrame(frameRight,100,200);
	 ButtonsFrame = new TGVerticalFrame(frameRight,100,200);
   // Create canvas widget
   fEcanvasID = new TRootEmbeddedCanvas("EcanvasID",frameLeft,100,100);
   frameLeft->AddFrame(fEcanvasID, new TGLayoutHints(kLHintsExpandX |
                   kLHintsExpandY, 10,10,10,1));

	 fEcanvasOD = new TRootEmbeddedCanvas("EcanvasOD",frameLeft,100,100);
	 frameLeft->AddFrame(fEcanvasOD, new TGLayoutHints(kLHintsExpandX |
									 kLHintsExpandY, 10,10,10,1));

   // Create a horizontal frame widget with buttons
   TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
	 TGTextButton *prev = new TGTextButton(hframe,"&Prev");
   prev->Connect("Clicked()","MyMainFrame",this,"Prev()");
   hframe->AddFrame(prev, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
   TGTextButton *next = new TGTextButton(hframe,"&Next");
   next->Connect("Clicked()","MyMainFrame",this,"Next()");
   hframe->AddFrame(next, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
   TGTextButton *exit = new TGTextButton(hframe,"&Exit",
                                "gApplication->Terminate(0)");
   hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));

	 TGRadioButton *radioButton[3];
	 TGButtonGroup *br = new TGButtonGroup(ButtonsFrame,"Show Detector",kVerticalFrame);
	 radioButton[0] = new TGRadioButton(br, new TGHotString("&ID"));
	 radioButton[1] = new TGRadioButton(br, new TGHotString("&OD"));
	 radioButton[2] = new TGRadioButton(br, new TGHotString("&ID and OD"));
	 br->Show();

	 radioButton[0]->Connect("Clicked()","MyMainFrame",this,"IDOnly()");
	 radioButton[1]->Connect("Clicked()","MyMainFrame",this,"ODOnly()");
	 radioButton[2]->Connect("Clicked()","MyMainFrame",this,"IDOD()");



	 strings[0] = "Event\t" + OutString(ev);
	 strings[1] = "Energy\t" + OutString(ene) + "MeV";
	 strings[2] = "vtx [cm]\t" + OutString(vtxx) + ", " + OutString(vtxy) + ", " + OutString(vtxz);
	 strings[3] = "dir [norm]\t" + OutString(dirx) + ", " + OutString(diry) + ", " + OutString(dirz);
	 strings[4] = "---------ID-------------";
	 strings[5] = "PMTs Hit Digi\t"+OutString(pmthitdigiid);
	 strings[6] = "PMTs Hit Raw\t"+OutString(pmthitrawid);
	 strings[7] = "Digi Hits\t"+OutString(hitdigiid);
	 strings[8] = "Raw Hits\t"+OutString(hitrawid);
	 strings[9] = "---------OD-------------";
	 strings[10] = "PMTs Hit Digi\t"+OutString(pmthitdigiod);
	 strings[11] = "PMTs Hit Raw\t"+OutString(pmthitrawod);
	 strings[12] = "Digi Hits\t"+OutString(hitdigiod);
	 strings[13] = "Raw Hits\t"+OutString(hitrawod);

	 for (int gh = 0; gh < 14; gh++){

		 text[gh] = new TGLabel(EventDetailsFrame, strings[gh].c_str());
		 EventDetailsFrame->AddFrame(text[gh], new TGLayoutHints(kLHintsCenterX,1,1,1,1));
	 }


	 //frameRight->AddFrame(EventDetailsFrame, new TGLayoutHints(kLHintsCenterX|kLHintsCenterY,0,0,0,0));
	 //frameRight->AddFrame(br, new TGLayoutHints(kLHintsLeft,5,5,3,4));
	 ButtonsFrame->AddFrame(br, new TGLayoutHints(kLHintsRight,0,0,0,0));
	 frameRight->AddFrame(EventDetailsFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));
	 frameRight->AddFrame(ButtonsFrame, new TGLayoutHints(kLHintsRight,0,0,0,0));

	 fMain->AddFrame(frameLeft, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,
																						 10,10,10,1));
	 fMain->AddFrame(frameRight, new TGLayoutHints(kLHintsCenterX| kLHintsBottom,
																						 10,10,10,1));
	 fMain->AddFrame(hframe, new TGLayoutHints(kLHintsBottom,
																						 2,2,2,2));


   // Set a name to the main frame
   fMain->SetWindowName("HK Event Display");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());
	 //frameLeft->Resize(frameLeft->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();

	 // Initialize other variables and plots
	 Vision();
	 Active();
	 UpdateText();
	 if (!odOn && idOn){IDOnly();radioButton[2]->SetEnabled(kFALSE);radioButton[1]->SetEnabled(kFALSE);radioButton[0]->SetState(kButtonDown);}
	 if (odOn && !idOn){ODOnly();radioButton[2]->SetEnabled(kFALSE);radioButton[0]->SetEnabled(kFALSE);radioButton[1]->SetState(kButtonDown);}
	 if(!odOn && !idOn){gApplication->Terminate(0);}
	 if(odOn && idOn){radioButton[2]->SetState(kButtonDown);}

}

void MyMainFrame::UpdateText(){

	strings[0] = "Event\t" + OutString(ev);
	strings[1] = "Energy\t" + OutString(ene) + "MeV";
	strings[2] = "vtx [cm]\t" + OutString(vtxx) + ", " + OutString(vtxy) + ", " + OutString(vtxz);
	strings[3] = "dir [norm]\t" + OutString(dirx) + ", " + OutString(diry) + ", " + OutString(dirz);
	strings[4] = "---------ID-------------";
	strings[5] = "PMTs Hit Digi\t"+OutString(pmthitdigiid);
	strings[6] = "PMTs Hit Raw\t"+OutString(pmthitrawid);
	strings[7] = "Digi Hits\t"+OutString(hitdigiid);
	strings[8] = "Raw Hits\t"+OutString(hitrawid);
	strings[9] = "---------OD-------------";
	strings[10] = "PMTs Hit Digi\t"+OutString(pmthitdigiod);
	strings[11] = "PMTs Hit Raw\t"+OutString(pmthitrawod);
	strings[12] = "Digi Hits\t"+OutString(hitdigiod);
	strings[13] = "Raw Hits\t"+OutString(hitrawod);

	for (int gh = 0; gh < 14; gh++){

		text[gh]->ChangeText(strings[gh].c_str());
	}
}


void MyMainFrame::IDOnly(){
	frameLeft->ShowFrame(fEcanvasID);
	frameLeft->HideFrame(fEcanvasOD);

}
void MyMainFrame::ODOnly(){
	frameLeft->ShowFrame(fEcanvasOD);
	frameLeft->HideFrame(fEcanvasID);

}
void MyMainFrame::IDOD(){
	frameLeft->ShowFrame(fEcanvasID);
	frameLeft->ShowFrame(fEcanvasOD);

}
void MyMainFrame::Prev() {
	if (ev <= 0 ) {std::cout << "This is the first event!" <<std::endl;}
	else{
		ev--;
		Active();
		UpdateText();
	}
}
void MyMainFrame::Next() {
	if (ev >= nEvent - 1 ) {std::cout << "This is the last event!" <<std::endl;}
	else{
		ev++;
		Active();
		UpdateText();
	}
}

MyMainFrame::~MyMainFrame() {
   // Clean up used widgets: frames, buttons, layout hints
   fMain->Cleanup();
   delete fMain;
}
void HKED(std::string file) {

   // Popup the GUI...
   new MyMainFrame(file, gClient->GetRoot(), 800, 600);
}
