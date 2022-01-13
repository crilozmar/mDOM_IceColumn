#include "mdomRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include <ctime>
#include <sys/time.h>

#include "mdomAnalysisManager.hh"
#include "mdomMinimization.hh"

extern std::vector<std::pair<G4int,G4int> > AcumulateHits;
extern MdomMinimization Minimization;

extern MdomAnalysisManager gAnalysisManager;
extern G4String gfilename;
extern G4int gsimevents;
extern G4int gneutroncapture;
extern G4int gSNGun;
extern G4String	gHittype;


mdomRunAction::mdomRunAction(){}
mdomRunAction::~mdomRunAction(){}

void mdomRunAction::BeginOfRunAction(const G4Run*)
{
	AcumulateHits.clear();
// 	Open output data file
    
//	With this you dont delete the file but write in next line (coment gAnalysisManager.WriteHeader() if you are using this one!)
	//gAnalysisManager.datafile.open(gfilename, std::ios::out| std::ios::app); 
//	This kill the current file and start from the beggining
//	gAnalysisManager.datafile.open(gfilename, std::ios::out| std::ios::trunc); 
    G4String infofilename;
    G4String filename_gamma2MeV;
    G4String filename_gamma8MeV;
    G4int lastdotposition;
    for (unsigned i=gfilename.length()-1; i>0; --i) {
        G4String character = gfilename.at(i);
        if (character == ".") {
            lastdotposition = i;
            break;
        }
    }
    for (unsigned i=0; i<lastdotposition; ++i) {
        G4String character = gfilename.at(i);
        infofilename.append(character);
        filename_gamma2MeV.append(character);
        filename_gamma8MeV.append(character);
        }
    infofilename.append("_info.txt");
    filename_gamma2MeV.append("_gamma2MeV.txt");
    filename_gamma8MeV.append("_gamma8MeV.txt");
    
    if (gSNGun != 0) {
        gAnalysisManager.maininfofile.open(infofilename, std::ios::out| std::ios::trunc); 
    }
    gAnalysisManager.datafile.open(gfilename, std::ios::out| std::ios::trunc); 

    if (gSNGun == 2) {
        if ((gneutroncapture == 1) || (gneutroncapture == 3)) {
            gAnalysisManager.datafile_gammas2MeV.open(filename_gamma2MeV, std::ios::out| std::ios::trunc); 
        } 
        if ((gneutroncapture == 2) || (gneutroncapture == 3)) {
            gAnalysisManager.datafile_gammas8MeV.open(filename_gamma8MeV, std::ios::out| std::ios::trunc); 
        }
    }
    gAnalysisManager.WriteHeader();
}

void mdomRunAction::EndOfRunAction(const G4Run*)
{
// Building likelihood function
  /*
	G4double TotalHits[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	G4int coordinate;

	for ( int j=0; j<(G4int)AcumulateHits.size(); j++ ) {
		coordinate = AcumulateHits[j].first;
		TotalHits[coordinate] = TotalHits[coordinate] + AcumulateHits[j].second;	  
	}
	
	for ( int j=0; j<24; j++ ) {
		TotalHits[j] = TotalHits[j]/gsimevents;
		gAnalysisManager.TotHits[j] = TotalHits[j];
	}

	Minimization.BuildLikelihood(TotalHits);
	
	// Write output file after all runs
	gAnalysisManager.Write();
	*/
        if (gHittype == "collective") {
            gAnalysisManager.WriteAccept();
            //gAnalysisManager.testwritter();
            gAnalysisManager.ResetEvent();
        }
	// Close output data file
	gAnalysisManager.datafile.close();
  
	//G4cout << "DEBUG: " << gAnalysisManager.hitStats.size() << G4endl;
  
}

