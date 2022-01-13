#include "OMSimDetectorConstruction.hh"
#include "mdomPhysicsList.hh"
#include "mdomPrimaryGeneratorAction.hh"
#include "mdomRunAction.hh"
#include "mdomEventAction.hh"
#include "mdomTrackingAction.hh"
#include "mdomSteppingAction.hh"
#include "mdomSteppingVerbose.hh"
#include "mdomAnalysisManager.hh"
#include "mdomMinimization.hh"
#include "mdomStackingAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "OMSimLogger.hh"

#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"	//xxx

#include "argtable2.h"
#include <ctime>
#include <sys/time.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

//#include "Minuit2/Minuit2Minimizer.h"
//#include "Math/Functor.h"
//#include "Math/Factory.h"

#include <cmath>	// for abs() of doubles
#include "G4SystemOfUnits.hh"	// if compiling with Geant4.10
#include "G4Navigator.hh"

unsigned int	stats_buffer_max_size = 10;	// how many hits to keep in memory before purging to file in EndOfEventAction
unsigned int	stats_buffer_last_purge_at = 0;	// at what hits count was the hits file last written to
// std::vector<G4int>	stats_PMT_hit;
// std::vector<G4int>	stats_OM_hit;

// double hc_eVnm = 1239.84193; // h*c in eV * nm
G4double	gworldsize;
G4double gRadius; // cylindrical world
G4double gHeight; // cylindrical world
G4int	gDepthpos;
G4int gneutroncapture;

G4ThreeVector gLEDpos;
G4int gEmitter_mdom;
const char*     gActivateLED;
G4double gInnercolumnradius;
G4double gInnercolumn_pos_x;// hicetube x pos relative to dom radius
G4double gInnercolumn_pos_y;// hicetube y pos relative to dom radius
//G4double gInnercolumn_av_costheta; //Our standard ice has 0.9
G4double gInnercolumn_b_inv;
G4String glightprofile_file;
G4double gZshift; //center of shooting module when LEDs are activated
G4double gCylHigh; //height of the cylindrical part of the mDOM
G4int gbulk_ice;

G4double	gscintYield;
G4double	gscintTimeConst;
G4double	gscintSpectrum;
G4double	gTemperature;
G4double	gRefCone_angle;

G4bool gharness_ropes;
G4double	gmdomseparation;
G4int	gn_mDOMs;
G4bool gpointingdownLED;

G4int	gsimevents;
G4int 	gReconstruction;
G4double	 LogL;
G4String	ggunfilename;
G4String	gnufluxname;
G4String	gnubarfluxname;
G4String	gSunspectrum;
G4double	gDistance;
G4double	gBeamDiam;
G4double	gtheta;
G4double	gphi;
G4double        gwavelen;
G4String	gfilename;
G4String	gEvent2Reconstruct;
G4String	gHittype;
G4String	gQEfile;
G4int		gPMT;
G4int		gEnvironment;
G4bool		gVisual;
G4bool		gInteractive;
G4bool		gHeader;
G4bool		gQE;
G4bool		gQEweigh;
G4bool		gSun_nue;
G4bool		gSun_numutau;
G4double	gSNmeanEnergy;
G4bool		gfixmeanenergy;
G4double 	gfixalpha;

G4int		gGlass;
G4int		gGel;
G4int		gConeMat;
G4int		gHolderColor;
G4int		gDOM;
G4int 		gSNGun;

G4String	gcode;
std::vector<std::tuple<G4int,G4int,G4int> > AcumulateHits;

std::stringstream command;


// G4String	greffilename;

G4bool		gKillAll;
G4long		current_event_id;

G4double	gposX, gposY, gposZ;
G4double	px,	py,	pz;
G4double 	gtMin, 	gfMin;

struct timeval	gTime_Run_Start;
struct timeval	gTime_Run_End;
long randseed;

MdomAnalysisManager gAnalysisManager;
MdomMinimization Minimization;
G4Navigator* aNavigator = new G4Navigator();


G4UImanager* UI;

// void clearstats() {
// 	stats_PMT_hit.clear();
// 	stats_OM_hit.clear();
// }

std::vector<G4String> explode(G4String s, char d) {
	std::vector<G4String> o;
	int i,j;
	i = s.find_first_of("#");
	if (i == 0) return o;
 	while (s.size() > 0) {
		i = s.find_first_of(d);
		j = s.find_last_of(d);
		o.push_back(s.substr(0, i));
		if (i == j) {
			o.push_back(s.substr(j+1));
			break;
		}
		s.erase(0,i+1);
 	}
	return o;// o beinhaltet s ohne d
}
std::vector<G4String> explode(char* cs, char d) {
	std::vector<G4String> o;
	G4String s = cs;
	return explode(s,d);
}
std::vector<double> readColumnDouble (G4String fn, int col) {
	std::vector<double>	values;
	unsigned int c;
	double	a;
	c = col;
	std::ifstream	infile;
	std::vector<G4String> n;
	char l[256];
	G4String l2;
	infile.open(fn);
	while (infile.good() && !infile.eof()) {
		infile.getline(l,255);
		l2 = l;
		n = explode(l2,'\t');
		if (n.size()>=c) {
			a = atof(n.at(c-1));
			values.push_back(a);
		}
	}
	infile.close();

	return values;//values enth??lt den c. Wert aus fn (aus jeder Spalte,welche  nach 255 zeichen oder durch \n beendet wird?)
}
/*
int myreadfiles(G4String filename) { 
  std::ifstream infile;   
  int num = 0; // num must start at 0
  infile.open(filename);// file containing numbers in 2 columns
  std::vector<G4double> column1;
  std::vector<G4double> column2;
     if(infile.fail()) // checks to see if file opended 
    { 
       G4cout << "error reading file :" << filename << G4endl; 
      return 1; // no point continuing if the file didn't open...
    } 
       while(!infile.eof()) // reads file to end of *file*, not line
      { 
         infile >> column1[num]; // read first column number
         infile >> column2[num]; // read second column number

         ++num; // go to the next number

         // you can also do it on the same line like this:
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
    for (unsigned int u = 0; u <column1.size(); u++) {
        G4cout << column1.at(u) << G4endl;
    }
  infile.close(); 

  return 0; // everything went right.
} */

void PointSource(double XGun, double YGun, double ZGun, double tGun, double fGun) {
  
	command.str("");
	command << "/gps/pos/centre "<< XGun <<" "<< YGun <<" "<< ZGun <<" mm";
	UI->ApplyCommand(command.str());
	
	G4double pi = CLHEP::pi;
	G4double pibar = pi/180.0;
	px = -sin(tGun*pibar)*cos(fGun*pibar);
	py = -sin(tGun*pibar)*sin(fGun*pibar);
	pz = -cos(tGun*pibar);
	
	command.str("");
	command << "/gps/direction " << px <<" "<< py <<" "<< pz;
	UI->ApplyCommand(command.str());
	
}
/*
G4double FunctionToReconstruction(const G4double *xx) {
  
	G4double fun;
	
	const G4double xMin = gposX;
	const G4double yMin = gposY;
	const G4double zMin = gposZ;
	const G4double tMin = xx[0];
	const G4double fMin = xx[1];
	
	gtMin = tMin;
	gfMin = fMin;
	
	G4cout << "****************************************************" << G4endl;
	G4cout << "****************************************************" << G4endl;
	G4cout << "****************************************************" << G4endl;
	G4cout << "****************************************************" << G4endl;
	G4cout << "x -> " << xMin << ",y -> " << yMin << ",z -> " << zMin << ",theta -> " << tMin << ",phi -> " << fMin << G4endl;
	G4cout << "****************************************************" << G4endl;
	G4cout << "****************************************************" << G4endl;
	G4cout << "****************************************************" << G4endl;
	G4cout << "****************************************************" << G4endl;

	  
	PointSource(xMin,yMin,zMin,tMin,fMin);
	
	if (gsimevents > 0) {
		command.str("");
		command << "/run/beamOn " << gsimevents;
		UI->ApplyCommand(command.str());
	}
	
	gAnalysisManager.Write();
	
	fun = Minimization.Minimize();
	G4cout <<  " !!!!!!!!!!!!!" <<fun << G4endl;

	return fun;
	
}
*/

void ReconstructionOfEvents() {
/*
	ROOT::Math::Minimizer* min = 
	      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
	
	ROOT::Math::Functor f(&FunctionToReconstruction,2);
   
	
	double maxworld = gworldsize/pow(3,1./2.)*1000; //*1000 if world size is define in metters!
	double variable[5] = {gposX,gposY,gposZ,gtheta,gphi};
	double step[5] = {10,10,10,60,60};
	//double lower[5] = {200.,200.,200.,0,-180};
	//double upper[5] = {maxworld,maxworld,maxworld, 180, 180};
 
	min->SetFunction(f);
  
	// Set the free variables to be minimized!
	
	min->SetVariable(0,"xMin",variable[0], step[0]);
	min->SetVariable(1,"yMin",variable[1], step[1]);
	min->SetVariable(2,"zMin",variable[2], step[2]);
	min->SetVariable(3,"tMin",variable[3], step[3]);
	min->SetVariable(4,"fMin",variable[4], step[4]);
	min->SetLimitedVariable(0 , "xMin" , variable[0] , step[0] , lower[0], upper[0]); 
	min->SetLimitedVariable(1 , "yMin" , variable[1] , step[1] , lower[1], upper[1]);
	min->SetLimitedVariable(2 , "zMin" , variable[2] , step[2] , lower[2], upper[2]);
	
	//min->SetLimitedVariable(0 , "tMin" , variable[3] , step[3] , lower[3], upper[3]);
	//min->SetLimitedVariable(1 , "fMin" , variable[4] , step[4] , lower[4], upper[4]);
	min->SetVariable(0,"tMin",variable[3], step[3]);
	min->SetVariable(1,"fMin",variable[4], step[4]);

	min->SetMaxFunctionCalls(gReconstruction);
	min->SetMaxIterations(10000);
	min->SetTolerance(0.00001);
 
	min->Minimize(); 
	const double *xs = min->X();
	G4cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
	    << FunctionToReconstruction(xs) << G4endl;  
 
	delete min;
	*/
}

void PointSourceGPS() {
	double cos_phi, sin_phi, cos_theta, sin_theta;
	double rho, posX,posY,posZ;
	
	cos_phi = cos(gphi*deg);
	sin_phi = sin(gphi*deg);
	cos_theta = cos(gtheta*deg);
	sin_theta = sin(gtheta*deg);
	
	rho = gDistance*sin_theta;
	posX = rho*cos_phi;
	posY = rho*sin_phi;
	posZ = gDistance*cos_theta;	
	
	G4cout << gDistance << "\t" << posX <<" "<< posY <<" "<< posZ << G4endl;
	
	command.str("");
	command << "/gps/position "<< posX <<" "<< posY <<" "<< posZ <<" m";
	UI->ApplyCommand(command.str());
	/*
	command.str("");
	command << "/gps/energy " << (1239.84193 / gwavelen) << " eV ";
	UI->ApplyCommand(command.str());
	*/
}


void LED_Definition(OMSimDetectorConstruction* detector)
{
	if (gDOM != 0) {
		G4String mssg = "LEDs are only defined for mDOM so far!";
    	critical(mssg);
		return;
	}
	std::vector<G4Transform3D> lLEDTransformers = detector->mOpticalModule->mLEDTransformers;
	std::vector<std::vector<G4double>> lLED_AngFromSphere = detector->mOpticalModule->mLED_AngFromSphere; //rho,theta,phi of the LEDs with respect to the center of the corresponding spherical part of the vessel
	G4RotationMatrix lShootingModuleOrientation = detector->mOpticalModule->mPlacedOrientations.at(gEmitter_mdom);
	G4ThreeVector lShootingModulePos = detector->mOpticalModule->mPlacedPositions.at(gEmitter_mdom);
	gZshift = lShootingModulePos.getZ();
	gCylHigh = detector->mOpticalModule->mCylHigh;

	G4int lNrTotalLED = detector->mOpticalModule->mNrTotalLED;
    int nbOfLEDs = strlen(gActivateLED); //number of activates/deactivates LEDs. For example, this can be 0011 to activate LEDs (3,4), even if lNrTotalLED is for example 10. 

	G4double DistLED = 152.8*mm; //Distance from center of sphere to LED light source?
    //G4double totalDOMsize = 199 * mm;

    double phi, theta, rho;
	double phi_led, theta_led,  rho_led;
    double posX, posY, posZ;
	double xPrime, yPrime, zPrime;
	double xPrime2, yPrime2, zPrime2;
    G4ThreeVector lLEDpos; //position of the LED with respect to the module
	G4Transform3D lLEDTrans; //position and rotation of the LED with respect to the world
	G4ThreeVector lLEDpos_fromworld; //position of the LED with respect to the world
	G4ThreeVector rot1vector;
	G4Transform3D rot1trans;
	G4ThreeVector rot2vector;
	G4Transform3D rot2trans;

    for (int i = 1; i < nbOfLEDs; i++)
    {
        static int iter = 0;
        if (*(gActivateLED + i) != '0')
        {
            command.str("");
            command << "/gps/source/add " << ++iter;
            UI->ApplyCommand(command.str());
        }
        
        if (i == (nbOfLEDs - 1) && *gActivateLED == '0')
        {
            command.str("");
            command << "/gps/source/delete " << iter;
            UI->ApplyCommand(command.str());
        }
    }

//     for (uint i = 0; i < v_LED_Positions.size(); i++)
    for (int i = 0; i < nbOfLEDs; i++)
    {
        static int iter = 0;
        if (*(gActivateLED + i) != '0')
        {
            command.str("");
            command << "/gps/source/set " << iter++;
            UI->ApplyCommand(command.str());
            
            phi_led = (lLED_AngFromSphere.at(i)).at(2);
            theta_led = (lLED_AngFromSphere.at(i)).at(1);
			rho_led = (lLED_AngFromSphere.at(i)).at(0);
            
            G4double gmDOMTiltingAngle_x = 0;
            G4double gmDOMTiltingAngle_y = 0;
            G4double xa = gmDOMTiltingAngle_x * deg; //NOTE: not used now
            G4double ya = gmDOMTiltingAngle_y * deg;//NOTE: not used now
            
			//get the positions like: position of the LED + position of the module with the corresponding rotation matrix
            posX = lLEDTransformers.at(i).getTranslation().getX();
			posY = lLEDTransformers.at(i).getTranslation().getY();
			posZ = lLEDTransformers.at(i).getTranslation().getZ();
			lLEDpos = G4ThreeVector(posX, posY, posZ);


			lLEDTrans = detector->mOpticalModule->GetNewPosition(lShootingModulePos, lShootingModuleOrientation, lLEDpos, G4RotationMatrix());
			lLEDpos_fromworld = lLEDTrans.getTranslation();
			posX = lLEDpos_fromworld.getX();
			posY = lLEDpos_fromworld.getY();
			posZ = lLEDpos_fromworld.getZ();

            command.str("");
            command << "/gps/pos/type Point";
            UI->ApplyCommand(command.str());
            
            if (gHittype == "collective") { //This means a direct light profile simulation will be done
                command.str("");
                command << "/gps/ang/type user"; //biast = theta
                UI->ApplyCommand(command.str());
                command.str("");
                command << "/gps/hist/type theta"; //biast = theta
                UI->ApplyCommand(command.str());

                std::vector<double> led_profile_x;
                std::vector<double> led_profile_y;
                led_profile_x = readColumnDouble(glightprofile_file, 1);
                led_profile_y = readColumnDouble(glightprofile_file, 2);
                for (unsigned int u = 0; u <led_profile_x.size(); u++) {
                    command.str("");
                    //G4cout << led_profile_x.at(u) << " " << led_profile_y.at(u) << G4endl;
                    command << "/gps/hist/point " << led_profile_x.at(u) << " " << led_profile_y.at(u); //biast = theta
                    UI->ApplyCommand(command.str());
                }
                command.str("");
                command << "/gps/ang/surfnorm" ;
                UI->ApplyCommand(command.str());
                
            } else {
                command.str("");
                //command << "/gps/ang/type iso";
                command << "/gps/ang/type cos";
                UI->ApplyCommand(command.str());
            }
            if (gpointingdownLED == true) {
                
                command.str("");
                command << "/gps/pos/rot1 " << 0. << " " << 1. << " " << 0.;
                UI->ApplyCommand(command.str());
                
                command.str("");
                command << "/gps/ang/rot1 " << 0. << " " << 1. << " " << 0.;
                UI->ApplyCommand(command.str());
                
                command.str("");
                command << "/gps/pos/rot2 " << -1. << " " << 0. << " " << 0.;
                UI->ApplyCommand(command.str());
                
                command.str("");
                command << "/gps/ang/rot2 " << -1. << " " << 0. << " " << 0.;
                UI->ApplyCommand(command.str());
                
                
            } else {
				//build angles rot1 and rot2. See https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html

				//Rot1
				xPrime=-sin(phi_led);
				yPrime=cos(phi_led);
				zPrime=0;
				//this would be the vector for rot1. However, same as before, we need to rotate this vector with the proper module rotation (no translation this time!)
				rot1vector = G4ThreeVector(xPrime, yPrime, zPrime);
				rot1trans = detector->mOpticalModule->GetNewPosition(G4ThreeVector(), lShootingModuleOrientation, rot1vector, G4RotationMatrix());
				rot1vector = rot1trans.getTranslation();
				xPrime = rot1vector.getX();
				yPrime = rot1vector.getY();
				zPrime = rot1vector.getZ();				
		
                command.str("");
                command << "/gps/pos/rot1 " << xPrime << " " << yPrime << " " << zPrime ;
                UI->ApplyCommand(command.str());
				command.str("");
                command << "/gps/ang/rot1 " << xPrime << " " << yPrime << " " << zPrime ;
                UI->ApplyCommand(command.str());
             
			 	//Rot2
				xPrime2=cos(theta_led)*cos(phi_led);
				yPrime2=cos(theta_led)*sin(phi_led);
				zPrime2=-sin(theta_led);
				//again this would be rot2, so we rotate
				rot2vector = G4ThreeVector(xPrime2, yPrime2, zPrime2);
				rot2trans = detector->mOpticalModule->GetNewPosition(G4ThreeVector(), lShootingModuleOrientation, rot2vector, G4RotationMatrix());
				rot2vector = rot2trans.getTranslation();
				xPrime2 = rot2vector.getX();
				yPrime2 = rot2vector.getY();
				zPrime2 = rot2vector.getZ();	
               
                command.str("");
                command << "/gps/pos/rot2 " << xPrime2 << " " << yPrime2 << " " <<  zPrime2;
				command.str("");
                command << "/gps/ang/rot2 " << xPrime2 << " " << yPrime2 << " " <<  zPrime2;
                UI->ApplyCommand(command.str());
                
                G4double MaxAngle = 89.; //when too close to 90, give photons that directly hit the structure and do not propagate... photons with theta=90 are anyway weighed very low
                command.str("");
                command << "/gps/ang/maxtheta " << MaxAngle << " deg";
                UI->ApplyCommand(command.str());
            }
            
            
            command.str("");
            command << "/gps/pos/centre " << posX << " " << posY << " " << posZ << " mm";

            UI->ApplyCommand(command.str());
            
            command.str("");
            command << "/gps/particle opticalphoton";
            //command << "/gps/particle geantino";
            UI->ApplyCommand(command.str());
            

            //command.str("");
            //command << "/gps/source/intensity 1";
            //UI->ApplyCommand(command.str());
            
            command.str("");
            command << "/gps/energy " << (1239.84193 / gwavelen) << " eV";
            UI->ApplyCommand(command.str());
            

            
            command.str("");
            command << "";
            UI->ApplyCommand(command.str());


        }
    }
}


int mdom() {
	struct timeval time_for_randy;
	gettimeofday(&time_for_randy, NULL);

	randseed = time_for_randy.tv_sec+4294*time_for_randy.tv_usec;
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(randseed,3));

	G4RunManager* runManager = new G4RunManager;


	OMSimDetectorConstruction* detector;
	detector = new OMSimDetectorConstruction();
	runManager->SetUserInitialization(detector);

	G4VUserPhysicsList* physics = new mdomPhysicsList;
	runManager->SetUserInitialization(physics);

        /*
	#ifdef G4VIS_USE
 		G4VisManager* visManager = new G4VisExecutive;
 		visManager->Initialize();
 		visManager->SetVerboseLevel(0);
	#endif
	*/
        auto visManager = new G4VisExecutive;
        visManager->Initialize();
        visManager->SetVerboseLevel(0);
        
	G4VUserPrimaryGeneratorAction* gen_action = new mdomPrimaryGeneratorAction();
	runManager->SetUserAction(gen_action);
    
    G4UserStackingAction* stacking_action = new mdomStackingAction();
    runManager->SetUserAction(stacking_action);

	G4UserRunAction* run_action = new mdomRunAction();
	runManager->SetUserAction(run_action);

	G4UserEventAction* event_action = new mdomEventAction();
	runManager->SetUserAction(event_action);

 	G4UserTrackingAction* tracking_action = new mdomTrackingAction();
 	runManager->SetUserAction(tracking_action);

	G4UserSteppingAction* stepping_action = new mdomSteppingAction();
	runManager->SetUserAction(stepping_action);

	runManager->Initialize();

	UI = G4UImanager::GetUIpointer();
	
	UI->ApplyCommand("/tracking/verbose 0");
	UI->ApplyCommand("/control/verbose 1");
	UI->ApplyCommand("/run/verbose 2");
	UI->ApplyCommand("/event/verbose 0");

	
	if (gReconstruction > 0) {
		ReconstructionOfEvents();
	} else { 
	  
	// Runnig GPS with imput parameters
	//	This is just for point source
		command.str("");
		command << "/selectGun "<< gSNGun;
		UI->ApplyCommand(command.str());
		
		if (gSNGun == 0){
			// setting up source:
			if (gActivateLED == "0000000000"){
				//PointSourceGPS();
				//Standard GPS
				command.str("");
				command << "/control/execute " << ggunfilename;
				UI->ApplyCommand(command.str());
			} else {
				LED_Definition(detector);
			}
		}
	// opening user interface prompt and visualization
		if (gInteractive){
			int argumc = 1;
			char* argumv[] = {"dummy", NULL};
			G4UIExecutive* UIEx = new G4UIExecutive(argumc, argumv);
			if (gVisual){
				UI->ApplyCommand("/control/execute ../aux/init_vis.mac"); //Warning!: VerboseLevels are redefined in init_vis.mac!!!
			}
			if (gsimevents > 0) {
			command.str("");
			command << "/run/beamOn " << gsimevents;
			UI->ApplyCommand(command.str());
			}
			UIEx->SessionStart();
			delete UIEx;
		}
		else{
			if (gsimevents > 0) {
			command.str("");
			command << "/run/beamOn " << gsimevents;
			UI->ApplyCommand(command.str());
			}
		}
	}

#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;
    if (gcode != "None") {
        G4cout << "Simulation finished : " << gcode << G4endl;
    }    
	return 0;
}


// G4int		gGlass;
// G4int		gGel;
// G4int		gConeMat;
// G4int		gHolderColor;
// G4int		gDOM;


int main(int argc,char *argv[])
{
	struct arg_dbl  *worldsize	= arg_dbl0("wW", "world","<n>","\t\tradius of world sphere in m");
	struct arg_dbl  *diameter	= arg_dbl0("dD", "diam","<n>","\t\tbeam diameter in mm");
	struct arg_dbl  *distance	= arg_dbl0("rR", "dist, rad","<n>","\t\temitter distance from surface of the mDOM, in m");
	struct arg_dbl  *xpos		= arg_dbl0("xX", "posx, posX","<n>","\t\t\tx of a puntual beam source in mm");
	struct arg_dbl  *ypos		= arg_dbl0("yY", "posy, posY","<n>","\t\t\ty of a puntual beam source in mm");
	struct arg_dbl  *zpos		= arg_dbl0("zZ", "posz, posZ","<n>","\t\t\tz of a puntual beam source in mm");
	struct arg_dbl  *theta		= arg_dbl0("tT", "theta","<n>","\t\ttheta (= zenith) in deg");
 	struct arg_dbl  *phi		= arg_dbl0("fF", "phi","<n>","\t\tphi (= azimuth) in deg");
	struct arg_dbl  *wavelen	= arg_dbl0("lL", "lambda","<n>","\t\twavelength of incoming light in nm");
	struct arg_int  *events		= arg_int0("nN", "numevents,nevents","<n>","\tnumber of photons emitted per angle");
	struct arg_file *gunfile	= arg_file0("gG","gun","<file.txt>","\t\tfile containing GPS parameters");
	struct arg_int  *pmt		= arg_int0("pP", "pmt,PMT","<n>","\t\tPMT type [12199S, etel, 12199e]"); 
        
        struct arg_str *activate_LED = arg_str0(NULL, "LED", "<e.g. 1001101010>", "\t\tactivating (1 - default) or deactivating (0) the LEDs (0 - 10)");
        struct arg_lit	*pointingdownLED		= arg_lit0(NULL,"pointingdownLED","\t\t Run the LEDs always pointing down, instead of being perpendicular to the surface (which is the default mode)");
	struct arg_int  *emitter_mdom		= arg_int0(NULL, "emitter_mdom","<n>","\t\tEmitter mDOM, default 0"); 
	struct arg_dbl  *innercolumn_radius		= arg_dbl0(NULL, "innercolumn_radius","<n>","\t\tRadius of inner ice column in cm, def. 7.5");
        struct arg_dbl  *innercolumn_pos_x		= arg_dbl0(NULL, "innercolumn_pos_x","<n>","\t\tx pos of the inner column relative to the modules, in cm. Def 0");
        struct arg_dbl  *innercolumn_pos_y		= arg_dbl0(NULL, "innercolumn_pos_y","<n>","\t\ty pos of the inner column relative to the modules, in cm. Def 0");
	//struct arg_dbl  *innercolumn_av_costheta		= arg_dbl0(NULL, "innercolumn_av_costheta","<n>","\t\tAv costheta for mie scattering in inner column"); 
	struct arg_dbl  *innercolumn_b_inv		= arg_dbl0(NULL, "innercolumn_b_inv","<n>","\t\tEffective scattering lenght for the inner column in cm. Def 3 cm"); 
    struct arg_str *lightprofile = arg_str0(NULL,"lightprofile","0 for direct simulation", "\t\tif != 0, it will directly simulates the chosen light experimental profile (from 4 to 11). Output will now be written after the whole run. [Def. 0]");
	struct arg_str  *extrap_mode		= arg_str0(NULL, "extrapmode","<n>","\t\t\tExtrapolation mode for lightprofile: 0(default): no extrapolation, 1:adds theta=90 with 0 intensity, 2:decrease intensity in 0.01 steps, 3:use the last 3 points to extrapolate until 0");
	struct arg_str *lplevel = arg_str0(NULL,"lplevel","<n>", "\t\tLevel of the light profile tu use [0:measured, 1:refraction correction, 2:whole correction (1+sin+binning)]. Default 2");
	struct arg_dbl  *mdomseparation		= arg_dbl0(NULL, "msep,mdomseparation","<n>","\t\t\tSeparation between mDOMs (center) in case that there is more than one in meters");
	struct arg_int  *n_mDOMs		= arg_int0(NULL, "nmdoms,n_mdoms","<n>","\t\tNumber of mDOMs in the simulation (0 = 1)"); 


        struct arg_int  *bulk_ice		= arg_int0(NULL, "bulk_ice","<n>","\t\tIceProperties in bulk ice (1) or not (0). If 0, --depthpos does nothing. Default should always be 1"); 
	struct arg_int  *glass		= arg_int0("uU", "glass","<n>","\t\t\tglass type [VITROVEX, Chiba, Kopp, myVitroVex, myChiba, WOMQuartz, fusedSilica]");
	struct arg_int	*gel 		= arg_int0("jJ", "gel", "<n>", "\t\t\tgel type [WACKER, Chiba, IceCube, Wacker_company]");
	struct arg_int	*conemat 	= arg_int0("kK", "conemat", "<n>", "\t\t\tcone material [V95, v98, aluminium, total98]");
	struct arg_int	*holdercol 	= arg_int0("cC", "holdercol", "<n>", "\t\t\tcone color [BLACK, white (Lambertian R = 98%)]");
	struct arg_int	*dom 		= arg_int0("mM", "om, dom", "<n>", "\t\t\tmodule type [MDOM, PDOM]");
	struct arg_dbl  *cone_ang   = arg_dbl0(NULL, "cone_ang","<n>","\t\t\topening semi-angle of cone; (45 deg)");	
	struct arg_int  *harness_ropes   = arg_int0(NULL, "harness_ropes","<n>","\t\t\tBuild complete harness and ropes if 1. Default = 1");	
        struct arg_str *code = arg_str0(NULL, "code", "None", "Whatever thing you want printout just before the GoodBye (None for no printing)");

	
	struct arg_dbl	*scintYield	= arg_dbl0(NULL, "scintYield", "<n>", "\t\tScintillation Yield of the glass (only Vitrovex). Default 57/MeV");
	struct arg_dbl	*scintTimeConst	= arg_dbl0(NULL, "scintTimeConst", "<n>", "\t\tScintillation's Time constant of the glass (only Vitrovex) in ns. Default 300000.");
	struct arg_dbl	*scintSpectrum	= arg_dbl0(NULL, "scintSpectrum", "<n>", "\t\tMove the scintillation's spectrum by # nm. Default 0 nm.");
	struct arg_dbl	*Temperature 	= arg_dbl0(NULL, "Temperature", "<n>", "\t\t Temperature for material property selection");
	
	struct arg_int  *environment= arg_int0("eE", "environment","<n>","\t\tmedium in which the setup is emmersed [AIR, ice, spice]");
	struct arg_file *outputfile	= arg_file0("oO","output","<file.txt>","\t\tfilename for hits data");
	struct arg_lit	*interactive= arg_lit0("iI","interact","\t\topens user interface after run");
	struct arg_lit	*visual		= arg_lit0("vV","visual","\t\tshows visualization of module after run (also calls interactive)");
	struct arg_lit	*nohead		= arg_lit0("qQ","nh, nohead","\t\tno header in outputfile");
	
	
	// My stuff
	struct arg_int  *depthpos	= arg_int0(NULL, "depthpos","<n>","\t\tDepth pos, check depth array in detectorconstruction to choose it propertly");
	struct arg_int *SN	= arg_int0(NULL,"SN","<n>","\t\t0=Heavy type II SN ls220, 1=Light type II SN ls220, 2=type IA SN DDT");
	struct arg_int	*SNGun 		= arg_int0(NULL, "SNgun", "<n>", "\t\tselect gun [GPS, SN ENEES, SN IBD [default without neutron capture]");
    struct arg_int  *neutroncapture      = arg_int0(NULL, "neutroncapture", "<n>", "\t\tif SNgun=IBD, choose whether simulate also the neutron capture or not [0=No neutron capture, 1=Neutron capture as 2 MeV gamma, 2=Neutron capture as 8 MeV gamma, 3=Neutron capture as 2 MeV and 8 MeV gammas (separated files)]");
	struct arg_lit	*Sun_e 		= arg_lit0(NULL, "Sun_e", "\t\tsimulate solar neutrinos -> nu_e interaction!");
	struct arg_lit	*Sun_mu 		= arg_lit0(NULL, "Sun_mu", "\t\tsimulate solar neutrinos -> nu_mu and nu_tau interaction!");
	struct arg_lit	*Sun_tau 		= arg_lit0(NULL, "Sun_tau", "\t\tsimulate solar neutrinos -> nu_mu and nu_tau interaction!");
	struct arg_dbl *SNmeanEnergy		= arg_dbl0(NULL, "SNmeanE","<n>","\t\tInstead of using SN model, use this mean energy");
	struct arg_dbl *alpha		= arg_dbl0(NULL, "alpha","<n>","\t\talpha parameter for the case of using --SNmeanEnergy");
	struct arg_int *reconstruction = arg_int0("aA","Reconstruction", "<n>", "\t\tnumber of calls for the reconstruction. By default 0 (no reconstruction)");
	struct arg_file *event2reconstruct = arg_file0("sS", "event2reconstruct", "<file.txt>","\t\tfile containing the hits per PMT to the event to reconstruct)");
	struct arg_file *QEfile = arg_file0(NULL,"QEfile"," <file.txt>","\t\tfile with the Quantum Efficiency");
	struct arg_lit *QE = arg_lit0(NULL,"QE","\t\tQuantum efficiency ON. Killing events");
	struct arg_lit *wQE = arg_lit0(NULL, "wQE", "\t\tWeigh of Quantum efficiency without killing the event");
	
	struct arg_lit	*help		= arg_lit0(NULL,"help","\t\tprint this help and exit");
	struct arg_end  *end		= arg_end(20);
	
	void* argtable[] = {worldsize,
                            diameter,
                            distance,
                            xpos, ypos, zpos,
                            theta, phi,
                            wavelen,
                            events,
                            gunfile,
                            pmt,
                            
                            activate_LED,
                            pointingdownLED,
                            emitter_mdom,
                            innercolumn_radius,
                            innercolumn_pos_x,
                            innercolumn_pos_y,
                            //innercolumn_av_costheta,
                            innercolumn_b_inv,
                            lightprofile,
							extrap_mode,
							lplevel,
                            
                            mdomseparation,
                            n_mDOMs,
                            
                            bulk_ice,
                            glass,
                            gel,
                            conemat,
                            holdercol,
                            dom,
                            cone_ang,
                            harness_ropes,
                            code,
                            
                            scintYield,
                            scintTimeConst,
                            scintSpectrum,
                            Temperature,
                            
                            environment,
                            outputfile,
                            interactive,
                            visual,
                            nohead,
                            depthpos,
                            SN, SNGun,neutroncapture,
                            Sun_e, Sun_mu, Sun_tau,
                            SNmeanEnergy,
                            alpha,
                            reconstruction,
                            event2reconstruct,
                            QEfile, QE, wQE,
                            help, end};

	const char* progname = "mdom";
	int nerrors;
	int exitcode=0;

	// verify the argtable[] entries were allocated sucessfully
	if (arg_nullcheck(argtable) != 0) {
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n",progname);
		exitcode=1;
		arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
		return exitcode;
	}
	// set any command line default values prior to parsing
	xpos->dval[0] = 0.0;
	ypos->dval[0] = 0.0;
	zpos->dval[0] = 300;
	worldsize->dval[0] = 10.0;	// world diameter in meters -- NOT USED!!!!!
	diameter->dval[0] = 420.0;	// 400 mm # for 14" sphere, 480 mm # for 17" sphere 
	distance->dval[0] = 1.0; // here value for mDOM scan with 2*mm safety margin
	theta->dval[0] = 90.0;
	phi->dval[0] = 0.0;
	wavelen->dval[0] = 405.0;	// [nm]
	events->ival[0] = 0;
	gunfile->filename[0] = "mdom.gps";
	pmt->ival[0] = 0;			// use new R12199 version as default
	
	//activate_LED->sval[0] = "11111111111111";
	activate_LED->sval[0] = "0000000000"; //default no LEDs
	emitter_mdom->ival[0] = 0;
	innercolumn_radius->dval[0] = 7.5; //cm
	innercolumn_pos_x->dval[0] = 0; //cm
	innercolumn_pos_y->dval[0] = 0; //cm
	//innercolumn_av_costheta->dval[0] = 0.95; //
	innercolumn_b_inv->dval[0] = 3; //cm
	bulk_ice->ival[0] = 1; //by default, scattering and absorption in bulk ice
	lightprofile->sval[0] = "0"; // def. 0, which means no direct simulation
	extrap_mode->sval[0] = "0"; //no extrapolation by default
	lplevel->sval[0] = "2"; //final level by default
        
	mdomseparation-> dval[0] = 2.7; //m
	n_mDOMs->ival[0] = 1;
	
	glass->ival[0] = 0;	// use VITROVEX as default
	gel->ival[0] = 1;	// Chiba gel as default	
	conemat->ival[0] = 0;	// use Alemco V95 as default
	holdercol->ival[0] = 0;	// use classic black holder as default
	dom->ival[0] = 0;	// use mDOM as default
	
	cone_ang->dval[0] = 51.0; // [degrees]	
	
	harness_ropes->ival[0] = 1;
    code->sval[0] = "None"; // core_process for example
	
	scintYield->dval[0] = 57.0;
	scintTimeConst->dval[0] = 300000.0;
	scintSpectrum->dval[0] = 0.0;
	Temperature->dval[0] = -35;
	
	//outputfile->filename[0] = "../output/mdom_testoutput_scan_angular.txt";
	SNmeanEnergy->dval[0] = 0.0; //<0.1 means do not take this parameter into account
	alpha->dval[0] = 3.0; //default alpha when using SNmeanEnergy
	// My stuff
	depthpos->ival[0] = 75; //Pos that have always been used of cleanest ice. Dirtiest ice is pos = 66
	environment->ival[0] = 2;	// use spice as default
	SN->ival[0] = 0; // Use heavy SN as default
	SNGun->ival[0]=0; //gps by default
	neutroncapture->ival[0]=0; //no neutron capture by default
	reconstruction->ival[0] = 0; //No reconstruction by default
	outputfile->filename[0] = "../ana/data.txt";
	event2reconstruct->filename[0] = "event2reconstruct.cfg"; //Just because we need a default file. This non-real event is all 0. 
	//The file that we use needs this structure, every row represents a PMT, and the number in the row is the number of hits in that PMT. It will be modified to take into account time. 
	QEfile->filename[0] = "QuantumEfficiency.cfg";
	
	/* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    G4String lightprofile_str = lightprofile->sval[0];
	G4String lplevel_str = lplevel->sval[0];
	G4String extrap_mode_str = extrap_mode->sval[0];
	//check if the levels and extrapolation modes are right
	std::string::size_type sz_lp;   // alias of size_t
	std::string::size_type sz_em;   // alias of size_t
	int extrap_mode_int = std::stoi (extrap_mode->sval[0],&sz_em);
	int lplevel_int = std::stoi (lplevel->sval[0],&sz_lp);
	if ((lplevel_int < 0) || (lplevel_int > 2)){
		G4cout << "ERROR: Wrong input level for light profile (only 0,1,2)" << G4endl;
		goto hell;
	}
	if ((extrap_mode_int < 0) || (extrap_mode_int > 3)){
		G4cout << "ERROR: Wrong input extrapolation mode for light profile (only 0,1,2)" << G4endl;
		goto hell;
	}

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
	{
        printf("\nGEANT4 simulation of the mDOM\n");
        printf("\nUsage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        arg_print_glossary(stdout,argtable,"  %-25s %s\n");
        printf("\n");
        exitcode=0;
        goto hell;
	}

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
	{
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto hell;
	}

    /* special case: uname with no command line options induces brief help */
    if (argc==1)
	{
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=0;
        goto hell;
	}

//	assign command-line arguments to variables:
	gworldsize = worldsize->dval[0];
	gBeamDiam = diameter->dval[0];
	//gDistance = (0.5 * 356.0 + 1.0)/1000.0 + distance->dval[0];
	gDistance = distance->dval[0];
	//gworldsize = gDistance + 1. ; //meters
	gtheta 	= theta->dval[0];
	gphi 	= phi->dval[0];
	gwavelen = wavelen->dval[0];
	gsimevents = events->ival[0];
	ggunfilename = gunfile->filename[0];
	gPMT = pmt->ival[0];
        
	if (strlen(activate_LED->sval[0]) > 10) {
		G4cout << "FATAL ERROR - INPUT PARAMETER: LEDs max number is so far 10!" << G4endl;
		goto hell;
	} else {
		gActivateLED = activate_LED->sval[0];
	}
	if (emitter_mdom->ival[0] >= n_mDOMs->ival[0]) {
		G4cout << "FATAL ERROR - INPUT PARAMETER: emitter mDOM number should be < nmdoms" << G4endl;
		goto hell;
	} else {
		gEmitter_mdom = emitter_mdom->ival[0];
	}
	
	if (pointingdownLED->count > 0) {
		gpointingdownLED = true;
	} else {
		gpointingdownLED = false;
	}
	gInnercolumnradius = innercolumn_radius->dval[0]*cm;
	gInnercolumn_pos_x = innercolumn_pos_x->dval[0]*cm;
	gInnercolumn_pos_y = innercolumn_pos_y->dval[0]*cm;
	/*if ((innercolumn_av_costheta->dval[0] > 1) || (innercolumn_av_costheta->dval[0] < 0)) {
		G4cout << "FATAL ERROR: innercolumn_av_costheta must be between 0 and 1!" << G4endl;
		goto hell;
	} else {
		gInnercolumn_av_costheta = innercolumn_av_costheta->dval[0];
	}*/
	gInnercolumn_b_inv = innercolumn_b_inv->dval[0] * cm;
	if ((bulk_ice->ival[0] != 0) && (bulk_ice->ival[0]) != 1) {
		G4cout << "FATAL ERROR - INPUT PARAMETER: bulk_ice must be between 0 and 1!" << G4endl;
		goto hell;
	} else {
		gbulk_ice = bulk_ice->ival[0];
	}

	gmdomseparation= mdomseparation->dval[0]*m;
	gn_mDOMs = n_mDOMs->ival[0];
	
	gGlass = glass->ival[0];
	gGel = gel->ival[0];
	gConeMat = conemat->ival[0];
	gHolderColor = holdercol->ival[0];
	gDOM = dom->ival[0];
	
	gRefCone_angle = cone_ang->dval[0];		
	
	gscintYield = scintYield->dval[0];
	gscintTimeConst = scintTimeConst->dval[0];
	gscintSpectrum = scintSpectrum->dval[0];
	gTemperature = Temperature->dval[0];
	
	gDepthpos = depthpos->ival[0];
	gEnvironment = environment->ival[0];
	gfilename = outputfile->filename[0];
	gposX = xpos->dval[0];
	gposY = ypos->dval[0];
	gposZ = zpos->dval[0];
	gSunspectrum = "solarneutrinospectrum.cfg";
	
	gSNmeanEnergy = SNmeanEnergy->dval[0];
	gfixalpha = alpha->dval[0];

	if (harness_ropes->ival[0]==1) gharness_ropes = true; else gharness_ropes = false;
    
    gcode = code->sval[0];
	
	if ((Sun_e->count>0) || (Sun_mu->count>0) || (Sun_tau->count>0) ){
		if ((Sun_e->count>0) && ((Sun_mu->count>0) || (Sun_tau->count>0))) {
			G4cout << "FATAL ERROR - INPUT PARAMETER: choose between nu_e or nu_mu/tau!!!" << G4endl;
			goto hell;
		}
		gSNGun = 3;
		G4cout << "*** Simulation of solar neutrinos ***" << G4endl;
		if (Sun_e->count>0) {
			gSun_nue =true;
			gSun_numutau =false;
		} else {  
			gSun_nue =false;
			gSun_numutau =true;
		}
	} else { //SNgun3 is to run solar neutrino interactions!
		gSNGun = SNGun->ival[0];
	}
		
	if ((neutroncapture->ival[0] < 0) || (neutroncapture->ival[0]>3)) {
        G4cout << "INPUT ERROR -> neutroncapture must be between 0 and 3!!!" <<G4endl;
        goto hell;
    } else {
        gneutroncapture = neutroncapture->ival[0];
    }
	gReconstruction = reconstruction->ival[0];
	gEvent2Reconstruct = event2reconstruct->filename[0];
	gQEfile = QEfile->filename[0];
	
	if (gSNmeanEnergy < 0.1) {
		gfixmeanenergy = false;
		if (SN->ival[0] == 0) {
			gnufluxname = "Flux_Nu_Heavy_ls220.cfg";
			gnubarfluxname = "Flux_Nubar_Heavy_ls220.cfg";
		} else if (SN->ival[0] == 1) {
			gnufluxname = "Flux_Nu_light_ls220.cfg";
			gnubarfluxname = "Flux_Nubar_light_ls220.cfg";
		} else if (SN->ival[0] == 2) {
			gnufluxname = "nu_DDT.cfg";
			gnubarfluxname = "nubar_DDT.cfg";
		} else if (SN->ival[0] == 3) {
			gnufluxname = "nu_GCD.cfg";
			gnubarfluxname = "nubar_GCD.cfg";
		} else {
			G4cout << "ERROR - INPUT PARAMETER!! Choose a valid SN model" << G4endl;
			goto hell;
		}
	} else {
		gfixmeanenergy = true;
		gnufluxname = "NoFile";
		gnubarfluxname = "NoFile";
	}
	
	if ((QE->count > 0) && (wQE->count == 0)) {
		gQE = true;
		gQEweigh = false;
	} else if ((QE->count == 0) && (wQE->count > 0)) {
		gQE = false;
		gQEweigh = true;
	} else if ((QE->count == 0) && (wQE->count == 0)) {
		gQE = false;
		gQEweigh = false;
	} else {
		G4cout << "FATAL ERROR - INPUT PARAMETER : can not choose QE killing event (--QE) and just weighing without killing (--wQE) at the same time" << G4endl;
		goto hell;
	}
	

	// Check if volume is big enough
	if (((gn_mDOMs % 2 == 0) && ((gmdomseparation*(gn_mDOMs/2-1./2.)) > (gworldsize*m+0.3*m))) || ((gn_mDOMs % 2 != 0) && ((gmdomseparation*(gn_mDOMs/2)) > (gworldsize*m+0.3*m)))) {
		G4cout << "FATAL ERROR - INPUT PARAMETER : This world is not big enough for those many mDOMs..." << G4endl;
		goto hell;
	}


	if (lightprofile_str == "0") {
		gHittype = "individual";
	} else {
		std::string::size_type sz;   // alias of size_t
		int lightprofile_int = std::stoi (lightprofile->sval[0],&sz);
		if ((lightprofile_int != 0) && ((lightprofile_int < 4) || (lightprofile_int > 13))) {
			G4cout << "FATAL ERROR - INPUT PARAMETER : --lightprofile only from 4 to 13 " << G4endl;
			goto hell;
			}
		gHittype = "collective";
		glightprofile_file = "../aux/processedlightspectrum_" + lightprofile_str + "_level_" + lplevel_str + "_ext" + extrap_mode_str + ".cfg";
		G4cout << "Using profile: " << glightprofile_file << G4endl;
	}
	
	if (interactive->count > 0) gInteractive = true; else gInteractive = false;
	if (visual->count > 0) {
		gVisual = true;
		gInteractive = true;
	}
	else {
		gVisual = false;
	}
	if (nohead->count > 0) gHeader = false; else gHeader = true;
	
	//	check params for sanity
	mdom();

hell:
    /* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
	//	return exitcode;

}


