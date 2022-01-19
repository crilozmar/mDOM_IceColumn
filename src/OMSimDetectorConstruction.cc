/** @file OMSimDetectorConstruction.cc
 *  @brief User defined detector.
 *
 * You should define and construct the detector here...this template is an example for a single arg module. 
 *
 *  @author 
 *  @date October 2021
 * 
 *  @version Geant4 10.7
 *  
 *  @todo 
 */

#include "OMSimDetectorConstruction.hh"
#include "abcDetectorComponent.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"

#include "OMSimInputData.hh"
#include "OMSimMDOM.hh"
#include "OMSimPDOM.hh"
#include "G4Navigator.hh"

extern G4double	gmdomseparation;
extern G4int	gn_mDOMs;
extern G4double gRadius; 
extern G4double gHeight; 
extern G4Navigator* aNavigator;
extern G4double gworldsize;
extern G4double gInnercolumn_pos_x;
extern G4double gInnercolumn_pos_y;
extern G4bool gharness_ropes;
//G4double gmDOMTiltingAngle_x = 0; // tilting angle of the mDOM w.r.t. the hole ice around x-axis
//G4double gmDOMTiltingAngle_y = 0; // tilting angle of the mDOM w.r.t. the hole ice around y-axis
extern G4double gInnercolumnradius;
//extern G4double gInnercolumn_av_costheta;


OMSimDetectorConstruction::OMSimDetectorConstruction()
    : mWorldSolid(0), mWorldLogical(0), mWorldPhysical(0)
{
}

OMSimDetectorConstruction::~OMSimDetectorConstruction()
{
}

/**
 * Construct the world volume
 */
void OMSimDetectorConstruction::ConstructWorld(){
    G4double innerRadiusOfTheTube = 0.*cm;
    G4double startAngleOfTheTube = 0.*deg;
    G4double spanningAngleOfTheTube = 360.*deg;
    gRadius = gworldsize*m;
    gHeight = gworldsize*m;
   
    mWorldSolid = new G4Tubs("World",
                 innerRadiusOfTheTube, 
                 gRadius,
                 gHeight,
                 startAngleOfTheTube, 
                 spanningAngleOfTheTube); 
  
    
    mWorldLogical = new G4LogicalVolume(mWorldSolid, mData->GetMaterial("argWorld"), "World_log", 0, 0, 0);
    //mWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4VisAttributes* World_vis= new G4VisAttributes(G4Colour(0.45,0.5,0.35,0.2));
    mWorldLogical->SetVisAttributes(World_vis);
    mWorldPhysical = new G4PVPlacement (0, G4ThreeVector(0.,0.,0.), mWorldLogical, "World_phys", 0, false, 0);
    aNavigator->SetWorldVolume(mWorldPhysical);

}

/**
 * Construct all solids needed for your study. Call OMSimOMConstruction for simulations with optical modules
 * and OMSimPMTConstruction for simulations with only PMTs.
 * @return World physical for the main
 */
G4VPhysicalVolume* OMSimDetectorConstruction::Construct() {
    bool check_overlap = true;

    mData = new OMSimInputData();
    mData->SearchFolders();
    
    ConstructWorld();
    //mData->SetValue("om_mDOM", "jNrPolarPMTs", 3);

    mOpticalModule = new mDOM(mData, gharness_ropes);
    
    //mOpticalModule->PlaceIt(G4ThreeVector(0,0,0), G4RotationMatrix(), mWorldLogical, "module_phys_0");

    bool lHoleIce = true;
    G4Tubs* HoleIceTube_inner_solid = new G4Tubs("HoleIceTube_inner_solid", 0, gInnercolumnradius, gworldsize*m, 0, 2 * CLHEP::pi);
    G4SubtractionSolid* HoleIceTube_inner_cut_solid;
    G4RotationMatrix lRot = G4RotationMatrix();

    

     if (gn_mDOMs <= 1) {
        mOpticalModule->PlaceIt(G4ThreeVector(0,0,0), lRot, mWorldLogical, "module_phys_0");
        G4Transform3D testrans = G4Transform3D(lRot, G4ThreeVector(0,0,0));
        if (lHoleIce) {
            HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(HoleIceTube_inner_solid, G4ThreeVector(-gInnercolumn_pos_x,-gInnercolumn_pos_y, 0), lRot, "HoleIceTube_inner_cut_solid_0");
        }
    } else {
        std::stringstream moduleconverter;
        std::stringstream cutglassconverter;
        G4bool firstcut = true;
        for (unsigned int k = 0; k < gn_mDOMs; k++){
            moduleconverter.str("");
            moduleconverter << "module_phys_" << k ;
            cutglassconverter.str("");
            cutglassconverter << "HoleIceTube_inner_cut_solid_" << k ; 
            G4cout << moduleconverter.str() << G4endl;
            G4double zpos;
            zpos = gmdomseparation*(gn_mDOMs/2.-k-1./2.);
            //G4cout << "Position z -> "<< zpos/m << G4endl;
            //G4cout << "Orientation -> "<< lRot.getPhi()/deg << ", "<< lRot.getTheta()/deg << ", "<<lRot.getPsi()/deg << ", "<< G4endl;
            //lRot.rotateY(23*deg);
            //lRot.rotateX(45*deg);
            //lRot.rotateZ(76*deg);
            mOpticalModule->PlaceIt(G4ThreeVector(0,0,zpos), lRot, mWorldLogical, moduleconverter.str());
            //G4Transform3D testrans = G4Transform3D(lRot, G4ThreeVector(0,0,zpos));
            //new G4PVPlacement(testrans, test, "test"+moduleconverter.str(), mWorldLogical, false, 0);
            if (lHoleIce) {
                if (firstcut == true) {
                    //HoleIceTube_inner_cut_solid = new G4SubtractionSolid(cutglassconverter.str(), HoleIceTube_inner_solid, lOmSolid, lRot, G4ThreeVector(-gInnercolumn_pos_x, -gInnercolumn_pos_y, zpos));
                    HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(HoleIceTube_inner_solid, G4ThreeVector(-gInnercolumn_pos_x,-gInnercolumn_pos_y, zpos), lRot, cutglassconverter.str());
                    firstcut = false;
                } else {
                    //HoleIceTube_inner_cut_solid = new G4SubtractionSolid(cutglassconverter.str(), HoleIceTube_inner_cut_solid, lOmSolid, lRot?, G4ThreeVector(-gInnercolumn_pos_x, -gInnercolumn_pos_y, zpos));
                    HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(HoleIceTube_inner_cut_solid, G4ThreeVector(-gInnercolumn_pos_x,-gInnercolumn_pos_y, zpos), lRot, cutglassconverter.str());
                }
            }
        }
    }
    

    if (lHoleIce) {
        G4LogicalVolume* HoleIceTube_inner_cut_logical = new G4LogicalVolume(HoleIceTube_inner_cut_solid,  mData->GetMaterial("Mat_BubColumn"), "HoleIceTube_inner_cut_logical");
        G4cout << gInnercolumn_pos_x/m << " , " << gInnercolumn_pos_y/m << G4endl;
        new G4PVPlacement(0, G4ThreeVector(gInnercolumn_pos_x, gInnercolumn_pos_y,0), HoleIceTube_inner_cut_logical, "HoleIceTube_inner_cut_physical", mWorldLogical, false, 0, check_overlap); 
        
        G4VisAttributes* HoleIceTube_inner_vis;
        HoleIceTube_inner_vis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.3, 0.4));
        HoleIceTube_inner_cut_logical->SetVisAttributes(HoleIceTube_inner_vis);   
    }
    return mWorldPhysical;
}
