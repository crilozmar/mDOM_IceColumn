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

extern G4double gdisplace_theta0;
extern G4double gdisplace_phi0;
extern G4double gdisplace_x0;
extern G4double gdisplace_y0;
extern G4double gdisplace_theta1;
extern G4double gdisplace_phi1;
extern G4double gdisplace_x1;
extern G4double gdisplace_y1;
extern G4double gtransitioncolumn;
extern G4int gtransitioncolumn_slices;


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
    G4Tubs* HoleIceTube_inner_solid;
    std::vector<G4VSolid*> columnsholder_solid;
    G4SubtractionSolid* HoleIceTube_inner_cut_solid;
    //std::vector<G4SubtractionSolid*> columnsholder_solid_cut;
    
    std::stringstream solidhole_converter;
    if (gtransitioncolumn != 0) {
        G4double start = gInnercolumnradius - gtransitioncolumn/2.;
        G4double end = gInnercolumnradius + gtransitioncolumn/2.;
        for (int k = 0; k < gtransitioncolumn_slices; k++){
            solidhole_converter.str("");
            solidhole_converter << "HoleIceTube_inner_solid_" << k;
            G4double rad = (end - start) / (gtransitioncolumn_slices-1) * k + start;
            //HoleIceTube_inner_solid = new G4Tubs(solidhole_converter.str(), 0, rad, gworldsize*m, 0, 2 * CLHEP::pi);
            columnsholder_solid.push_back( new G4Tubs(solidhole_converter.str(), 0, rad, gworldsize*m, 0, 2 * CLHEP::pi) );
        }
    } else {
        HoleIceTube_inner_solid = new G4Tubs("HoleIceTube_inner_solid", 0, gInnercolumnradius, gworldsize*m, 0, 2 * CLHEP::pi);
    }
    
    G4RotationMatrix lRot;
    G4double xpos;
    G4double ypos;
    

    if (gn_mDOMs <= 1) {
        lRot = G4RotationMatrix();
        mOpticalModule->PlaceIt(G4ThreeVector(0,0,0), lRot, mWorldLogical, "module_phys_0");
        G4Transform3D testrans = G4Transform3D(lRot, G4ThreeVector(0,0,0));
        if (lHoleIce) {
            if (gtransitioncolumn != 0) {
                for (int k = 0; k < gtransitioncolumn_slices; k++){
                    solidhole_converter.str("");
                    solidhole_converter << "HoleIceTube_inner_cut_solid_" << k;
                    //HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(columnsholder_solid.at(k), G4ThreeVector(-gInnercolumn_pos_x,-gInnercolumn_pos_y, 0), lRot, solidhole_converter.str());
                    columnsholder_solid.at(k) = mOpticalModule->SubstractToVolume(columnsholder_solid.at(k), G4ThreeVector(-gInnercolumn_pos_x,-gInnercolumn_pos_y, 0), lRot, solidhole_converter.str()) ;
                }
            } else {
                HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(HoleIceTube_inner_solid, G4ThreeVector(-gInnercolumn_pos_x,-gInnercolumn_pos_y, 0), lRot, "HoleIceTube_inner_cut_solid_0");
            }
        }
    } else {
        std::stringstream moduleconverter;
        std::stringstream cutglassconverter;
        G4bool firstcut = true;
        for (unsigned int k = 0; k < gn_mDOMs; k++){
            lRot = G4RotationMatrix();
            moduleconverter.str("");
            moduleconverter << "module_phys_" << k ;
            G4cout << moduleconverter.str() << G4endl;
            G4double zpos;
            zpos = gmdomseparation*(gn_mDOMs/2.-k-1./2.);
            //G4cout << "Position z -> "<< zpos/m << G4endl;
            //G4cout << "Orientation -> "<< lRot.getPhi()/deg << ", "<< lRot.getTheta()/deg << ", "<<lRot.getPsi()/deg << ", "<< G4endl;
            //lRot.rotateY(23*deg);
            //lRot.rotateX(45*deg);
            //lRot.rotateZ(76*deg);
            if (k == 0) {
                xpos = gdisplace_x0;
                ypos = gdisplace_y0;
                if (gdisplace_theta0 < 0.*deg) {
                    lRot.setTheta(gdisplace_theta0);
                    lRot.setPhi(180*deg + gdisplace_phi0);
                } else {
                    lRot.setTheta(gdisplace_theta0);
                    lRot.setPhi(gdisplace_phi0);
                }
            } else if (k == 1) {
                xpos = gdisplace_x1;
                ypos = gdisplace_y1;
                if (gdisplace_theta1 < 0.*deg) {
                    lRot.setTheta(gdisplace_theta1);
                    lRot.setPhi(180*deg + gdisplace_phi1);
                } else {
                    lRot.setTheta(gdisplace_theta1);
                    lRot.setPhi(gdisplace_phi1);
                }
            } else {
                xpos = 0*cm;    
                ypos = 0*cm;    
            }
            G4cout << "teeeeeeeest" << G4endl;
            mOpticalModule->PlaceIt(G4ThreeVector(xpos,ypos,zpos), lRot, mWorldLogical, moduleconverter.str());
            //G4Transform3D testrans = G4Transform3D(lRot, G4ThreeVector(0,0,zpos));
            //new G4PVPlacement(testrans, test, "test"+moduleconverter.str(), mWorldLogical, false, 0);
            if (lHoleIce) {
                if (gtransitioncolumn != 0) {
                                G4cout << "teeeeeeeest" << G4endl;
                    for (int u = 0; u < gtransitioncolumn_slices; u++){
                        cutglassconverter.str("");
                        cutglassconverter << "HoleIceTube_inner_cut_solid_" << u <<"_" << k;
                        G4cout << cutglassconverter.str() << G4endl;
                            HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(columnsholder_solid.at(u), G4ThreeVector(xpos-gInnercolumn_pos_x,ypos-gInnercolumn_pos_y, zpos), lRot, cutglassconverter.str());
                            columnsholder_solid.at(u) = HoleIceTube_inner_cut_solid;
                    }
                } else {
                    if (firstcut == true) {
                        //HoleIceTube_inner_cut_solid = new G4SubtractionSolid(cutglassconverter.str(), HoleIceTube_inner_solid, lOmSolid, lRot, G4ThreeVector(-gInnercolumn_pos_x, -gInnercolumn_pos_y, zpos));
                        HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(HoleIceTube_inner_solid, G4ThreeVector(xpos-gInnercolumn_pos_x,ypos-gInnercolumn_pos_y, zpos), lRot, cutglassconverter.str());
                        firstcut = false;
                    } else {
                        //HoleIceTube_inner_cut_solid = new G4SubtractionSolid(cutglassconverter.str(), HoleIceTube_inner_cut_solid, lOmSolid, lRot?, G4ThreeVector(-gInnercolumn_pos_x, -gInnercolumn_pos_y, zpos));
                        HoleIceTube_inner_cut_solid = mOpticalModule->SubstractToVolume(HoleIceTube_inner_cut_solid, G4ThreeVector(xpos-gInnercolumn_pos_x,ypos-gInnercolumn_pos_y, zpos), lRot, cutglassconverter.str());
                    }
                }
            }
        }
    }
    /*
    G4Material* testeando = mData->GetMaterial("Mat_BubColumn_0");
    G4cout << "aaaaaaaaaaaaaaaaaaaaaaaaAA" << G4endl;
    G4cout << testeando->GetName() << G4endl;
    G4MaterialPropertiesTable * testable = testeando->GetMaterialPropertiesTable();
    testable->DumpTable();
    
    G4Material* testeando2 = mData->GetMaterial("Mat_BubColumn_1");
    G4cout << "aaaaaaaaaaaaaaaaaaaaaaaaAA" << G4endl;
    G4cout << testeando2->GetName() << G4endl;
    G4MaterialPropertiesTable * testable2 = testeando2->GetMaterialPropertiesTable();
    testable2->DumpTable();
    */
    std::vector<G4LogicalVolume*> columnsholder_logicals;
    G4VisAttributes* HoleIceTube_inner_vis;
    if (lHoleIce) {                
        if (gtransitioncolumn != 0) {
            std::stringstream matname;
            std::stringstream holelogicalsname;
            std::stringstream holephysicalsname;
            for (int u = 0; u < gtransitioncolumn_slices; u++){ //logicals
                matname.str("");
                matname << "Mat_BubColumn_" << u;
                holelogicalsname.str("");
                holelogicalsname << "HoleIceTube_inner_cut_logical_" << u;
                columnsholder_logicals.push_back( new G4LogicalVolume(columnsholder_solid.at(u),  mData->GetMaterial(matname.str()), holelogicalsname.str()) );
                HoleIceTube_inner_vis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.3, 0.4));
                columnsholder_logicals.at(u)->SetVisAttributes(HoleIceTube_inner_vis);   
            }
            for (int u = gtransitioncolumn_slices-1; u >= 0; u--){ //placement in reverse
                holephysicalsname.str("");
                holephysicalsname << "HoleIceTube_inner_cut_physical_" << u;
                if (u == gtransitioncolumn_slices-1) {
                    new G4PVPlacement(0, G4ThreeVector(gInnercolumn_pos_x, gInnercolumn_pos_y,0), columnsholder_logicals.at(u), holephysicalsname.str(), mWorldLogical, false, 0, check_overlap); 
                } else {
                    new G4PVPlacement(0, G4ThreeVector(gInnercolumn_pos_x, gInnercolumn_pos_y,0), columnsholder_logicals.at(u), holephysicalsname.str(), columnsholder_logicals.at(u+1), false, 0, check_overlap); 
                }
            }
        } else {
            G4LogicalVolume* HoleIceTube_inner_cut_logical = new G4LogicalVolume(HoleIceTube_inner_cut_solid,  mData->GetMaterial("Mat_BubColumn"), "HoleIceTube_inner_cut_logical");
            //G4cout << gInnercolumn_pos_x/m << " , " << gInnercolumn_pos_y/m << G4endl;
            new G4PVPlacement(0, G4ThreeVector(gInnercolumn_pos_x, gInnercolumn_pos_y,0), HoleIceTube_inner_cut_logical, "HoleIceTube_inner_cut_physical", mWorldLogical, false, 0, check_overlap); 
            
            HoleIceTube_inner_vis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.3, 0.4));
            HoleIceTube_inner_cut_logical->SetVisAttributes(HoleIceTube_inner_vis);   
        }
    }
    /*
    for (int u = 0; u < gtransitioncolumn_slices; u++){ 
        G4cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << G4endl;
        G4cout << columnsholder_logicals.at(u)->GetName() << G4endl;
    }*/
    return mWorldPhysical;
}
