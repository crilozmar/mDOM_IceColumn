#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "OMSimDataFileTypes.hh"
#include "OMSimInputData.hh"
#include "OMSimLogger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include <G4UnitsTable.hh>
#include <dirent.h>
#include <cmath>

namespace pt = boost::property_tree;
extern G4double gInnercolumn_b_inv;
extern G4double gtransitioncolumn;
extern G4int gtransitioncolumn_slices;
extern G4int gtransitioncolumn_type;
extern G4double gInnercolumnradius;
extern G4double gMieVar;
extern G4double gAbsVar;

/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                                  Base Abstract Classes
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
/**
 * Constructor. Saves file name to member.
 * @param pFileName
 */
abcDataFile::abcDataFile(G4String pFileName)
{   mFileData = new ParameterTable();
    mFileName = pFileName;
}

template <typename T>
/**
 * Transforms the values inside a pTree-array to a vector. The values can be also transformed to a G4double.
 * @param pVector  vector where the (transformed) values are saved
 * @param pTree    pTree containing json data
 * @param pKey     json attribute label where values are found
 * @param pScaling Values of array are multiplied by this factor. You can set the fisical unit with this.
 * @param pInverse In case you need the inverse of a value x, 1/x is appended (e.g. transforming from nm to eV)
 */
void abcDataFile::ParseToVector(std::vector<T> &pVector, pt::ptree pTree, std::basic_string<char> pKey, G4double pScaling, bool pInverse)
{
    for (pt::ptree::value_type &ridx : pTree.get_child(pKey))
    { //get array from element with key "pKey" of the json
        if (pInverse)
        { // if we need 1/x
            pVector.push_back(pScaling / ridx.second.get_value<T>());
        }
        else
        { // otherwise we only by scaling factor
            pVector.push_back(ridx.second.get_value<T>() * pScaling);
        }
    }
}

/**
 * Defines new material from data in json-file. 
 */
void abcMaterialData::CreateMaterial()
{

    pt::read_json(mFileName, mJsonTree); //read json file into mJsonTree

    mFileData->AppendParameterTable(mFileName);
    mMPT = new G4MaterialPropertiesTable();
    mMatDatBase = G4NistManager::Instance();

    mObjectName = mJsonTree.get<G4String>("jName");
    const G4String lDataType = mJsonTree.get<G4String>("jDataType");

    const G4double lDensity = mFileData->GetValue(mObjectName, "jDensity");

    const G4String lState_str = mJsonTree.get<G4String>("jState");
    const G4State lState = GetState(lState_str);

    //Defining the material with its density, number of components, state and name
    mMaterial = new G4Material(mObjectName, lDensity, mJsonTree.get_child("jComponents").size(), lState);

    //Construct material with fractional components (isotopes or G4-Materials)
    for (pt::ptree::value_type &key : mJsonTree.get_child("jComponents"))
    {
        std::string componentName = key.first;
        double componentFraction = key.second.get_value<double>();
        mMaterial->AddMaterial(mMatDatBase->FindOrBuildMaterial(componentName), componentFraction);
    }
    G4String mssg = "New Material defined: " + mMaterial->GetName();
    info(mssg);
}
/**
 * Extracts absorption length and adds it to the material property table 
 */
void abcMaterialData::ExtractAbsorptionLength()
{
    std::vector<G4double> lAbsLength;
    std::vector<G4double> lAbsLengthEnergy;
    ParseToVector(lAbsLength, mJsonTree, "jAbsLength", 1 * mm, false);
    ParseToVector(lAbsLengthEnergy, mJsonTree, "jAbsLengthWavelength", mHC_eVnm, true);
    mMPT->AddProperty("ABSLENGTH", &lAbsLengthEnergy[0], &lAbsLength[0], static_cast<int>(lAbsLength.size()));
}
/**
 * Extracts refraction index and adds it to the material property table 
 */
void abcMaterialData::ExtractRefractionIndex()
{
    std::vector<G4double> lRefractionIndex;
    std::vector<G4double> lRefractionIndexEnergy;
    ParseToVector(lRefractionIndex, mJsonTree, "jRefractiveIdx", 1., false);
    ParseToVector(lRefractionIndexEnergy, mJsonTree, "jRefractiveIdxWavelength", mHC_eVnm, true);
    mMPT->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
}
/**
 * State in string to G4State
 * @param  G4String
 * @return G4State
 */
G4State abcMaterialData::GetState(G4String pState_str)
{
    G4State lState;
    if (pState_str == "kStateSolid")
        lState = kStateSolid;
    else if (pState_str == "kStateLiquid")
        lState = kStateLiquid;
    else if (pState_str == "kStateGas")
        lState = kStateGas;
    else
        lState = kStateUndefined;
    return lState;
}
/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *                                     Derived Classes
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
/**
 * Extracts and creates material for material with refraction index and absorption length defined.
 */
void RefractionAndAbsorption::ExtractInformation()
{
    CreateMaterial();
    ExtractAbsorptionLength();
    ExtractRefractionIndex();
    mMaterial->SetMaterialPropertiesTable(mMPT);
}
/**
 * Extracts and creates material for material with refraction index defined.
 */
void RefractionOnly::ExtractInformation()
{
    CreateMaterial();
    ExtractRefractionIndex();
    mMaterial->SetMaterialPropertiesTable(mMPT);
}
/**
 * Extracts and creates material without optical properties.
 */
void NoOptics::ExtractInformation()
{
    CreateMaterial();
}
/**
 * Extracts and creates ice with optical properties from IceCube.
 */
void IceCubeIce::ExtractInformation()
{
    CreateMaterial(); //creates IceCubeICE
    G4Material *lIceMie = new G4Material("IceCubeICE_SPICE", mFileData->GetValue(mObjectName,"jDensity"), mMatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid); //create IceCubeICE_SPICE
    G4Material *lBubleColumnMie = new G4Material("Mat_BubColumn", mFileData->GetValue(mObjectName,"jDensity"), mMatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid); //create IceCubeICE_SPICE
    std::vector<G4Material*> lBubleColumnMie_slices;

    std::vector<G4double> lMieScatteringLength;
    std::vector<G4double> lMieScatteringLength_BubleColumn;
    std::vector<G4double> this_lMieScatteringLength_BubleColumn;
    std::vector<std::vector<G4double>> lMieScatteringLength_BubleColumn_slices;
    G4double this_Innercolumn_b_inv;
    std::vector<G4double> lWavelength;
    std::vector<G4double> lRefractionIndex;
    std::vector<G4double> lRefractionIndexEnergy;
    std::vector<G4double> lAbsLength;
    ParseToVector(lRefractionIndexEnergy, mJsonTree, "jWavelength_spice", mHC_eVnm, true);
    ParseToVector(lWavelength, mJsonTree, "jWavelength_spice", 1 * nm, false);
    ParseToVector(mSpice_be400inv, mJsonTree, "jbe400inv_spice", 1 * m, false);
    ParseToVector(mSpice_a400inv, mJsonTree, "ja400inv_spice", 1 * m, false);
    ParseToVector(mSpice_Depth, mJsonTree, "jDepth_spice", 1 * m, false);
    
    //loop to include variations on the effective scattering and absorption lenghts
    for (int u = 0; u < static_cast<int>(lRefractionIndexEnergy.size()); u++){
        mSpice_be400inv[u] = mSpice_be400inv.at(u) * (1.+ gMieVar/100.);
        mSpice_a400inv[u] = mSpice_a400inv.at(u) * (1. + gAbsVar/100.);
    }
    for (int u = 0; u < static_cast<int>(lRefractionIndexEnergy.size()); u++)
    {
        lRefractionIndex.push_back(Spice_Refraction(lWavelength.at(u)));
        lAbsLength.push_back(Spice_Absorption(lWavelength.at(u)));
        lMieScatteringLength.push_back(Mie_Scattering(lWavelength.at(u)));
        lMieScatteringLength_BubleColumn.push_back(gInnercolumn_b_inv);
    }
    if (gtransitioncolumn != 0) {
        for (int k = 0; k < static_cast<int>(gtransitioncolumn_slices); k++) {
            this_lMieScatteringLength_BubleColumn.clear();
            for (int u = 0; u < static_cast<int>(lRefractionIndexEnergy.size()); u++) {
                //G4cout << k << G4endl;
                if (gtransitioncolumn_type == 0) { //lineal
                    this_Innercolumn_b_inv = (Mie_Scattering(lWavelength.at(u)) - gInnercolumn_b_inv)/(gtransitioncolumn_slices-1) * k + gInnercolumn_b_inv; // index*difference+min
                } else { //gaussian
                    G4double A = Mie_Scattering(lWavelength.at(u));
                    G4double sigma = gtransitioncolumn/10.; //Dont forget this!!!! The larger gtransitioncolumn is used so that the geometry within the slices includes the whole transition
                    G4double mu = gInnercolumnradius; //so it is more or less centered
                    G4double y0 = gInnercolumn_b_inv;
                    //now assign, from a equally distribution in r, each slice k to one value of the gaussian
                    G4double start = gInnercolumnradius-gtransitioncolumn/2.;
                    G4double end = gInnercolumnradius+gtransitioncolumn/2.;
                    G4double rad = (end - start) / (gtransitioncolumn_slices-1) * k + start;
                    //now calculate the gaussian, where x=rad
                    //this_Innercolumn_b_inv = (A-y0) * exp(-pow((rad - mu)/sigma,2)/2.0) + y0;
                    this_Innercolumn_b_inv = 1./2.*(1.+erf((rad-mu)/(sigma*sqrt(2))))*A+y0;
                }
                this_lMieScatteringLength_BubleColumn.push_back(this_Innercolumn_b_inv);
                /*
                if (u == 30) {
                    G4cout << Mie_Scattering(lWavelength.at(u)) << ", " <<  gInnercolumn_b_inv << ", " <<  gtransitioncolumn_slices << ", " <<  k << ", " << gInnercolumn_b_inv << " : " << this_Innercolumn_b_inv << G4endl;
                }*/
            }
            lMieScatteringLength_BubleColumn_slices.push_back(this_lMieScatteringLength_BubleColumn);
        }
    }
    /*
    if (gtransitioncolumn != 0) {
                    int u = 30;
            G4cout << "Wavelenght : " << lWavelength.at(u)/nm << G4endl;
            G4cout << "be out of the column: " << lMieScatteringLength.at(u) << G4endl;
        for (int k = 0; k < static_cast<int>(lMieScatteringLength_BubleColumn_slices.size()); k++)
        {

            //G4cout << lMieScatteringLength_BubleColumn_slices.size() << " " << lMieScatteringLength_BubleColumn_slices[u].size() << G4endl;
            //for (int k = 0; k < static_cast<int>(lMieScatteringLength_BubleColumn_slices[u].size()); k++) {
                G4cout << k << G4endl;
                G4cout << "be array: " << lMieScatteringLength_BubleColumn_slices.at(k).at(u) << G4endl;
            //}
        }
    }
    */
    
    //give refractive index to IceCubeICE. This is used also for IceCubeICE_SPICE
    mMPT->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
    mMaterial->SetMaterialPropertiesTable(mMPT);
    
    //give properties to IceCubeICE_SPICE
    mMPT->AddProperty("ABSLENGTH", &lRefractionIndexEnergy[0], &lAbsLength[0], static_cast<int>(lAbsLength.size()));
    mMPT->AddProperty("MIEHG", &lRefractionIndexEnergy[0], &lMieScatteringLength[0], static_cast<int>(lRefractionIndex.size()))->SetSpline(true);
    mMPT->AddConstProperty("MIEHG_FORWARD", mMIE_spice_const[0]);
    mMPT->AddConstProperty("MIEHG_BACKWARD", mMIE_spice_const[1]);
    mMPT->AddConstProperty("MIEHG_FORWARD_RATIO", mMIE_spice_const[2]);
    lIceMie->SetMaterialPropertiesTable(mMPT);
    G4String mssg = "Ice properties at depth " + std::to_string(mSpice_Depth[mSpiceDepth_pos] / m) + " m.";
    notice(mssg);
    //now give the properties to the bubble column, which are basically the same ones but with the chosen scattering lenght
    mMPT_holeice = new G4MaterialPropertiesTable();
    mMPT_holeice->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
    mMPT_holeice->AddProperty("ABSLENGTH", &lRefractionIndexEnergy[0], &lAbsLength[0], static_cast<int>(lAbsLength.size()));
    mMPT_holeice->AddProperty("MIEHG", &lRefractionIndexEnergy[0], &lMieScatteringLength_BubleColumn[0], static_cast<int>(lRefractionIndex.size()))->SetSpline(true);
    mMPT_holeice->AddConstProperty("MIEHG_FORWARD", mMIE_spice_const[0]);
    mMPT_holeice->AddConstProperty("MIEHG_BACKWARD", mMIE_spice_const[1]);
    mMPT_holeice->AddConstProperty("MIEHG_FORWARD_RATIO", mMIE_spice_const[2]);
    lBubleColumnMie->SetMaterialPropertiesTable(mMPT_holeice);
    
    std::stringstream matname;
    G4MaterialPropertiesTable* MPT_holeice_temp;
    G4Material* lBubleColumnMie_temp;
    if (gtransitioncolumn != 0) {
        for (int k = 0; k < static_cast<int>(gtransitioncolumn_slices); k++) {
            matname.str("");
            matname << "Mat_BubColumn_" << k;
            this_lMieScatteringLength_BubleColumn.clear();
            this_lMieScatteringLength_BubleColumn = lMieScatteringLength_BubleColumn_slices.at(k);
            lBubleColumnMie_temp = new G4Material(matname.str(), mFileData->GetValue(mObjectName,"jDensity"), mMatDatBase->FindOrBuildMaterial("G4_WATER"), kStateSolid); //create IceCubeICE_SPICE
            
            MPT_holeice_temp = new G4MaterialPropertiesTable();
            MPT_holeice_temp->AddProperty("RINDEX", &lRefractionIndexEnergy[0], &lRefractionIndex[0], static_cast<int>(lRefractionIndex.size()));
            MPT_holeice_temp->AddProperty("ABSLENGTH", &lRefractionIndexEnergy[0], &lAbsLength[0], static_cast<int>(lAbsLength.size()));
            MPT_holeice_temp->AddProperty("MIEHG", &lRefractionIndexEnergy[0], &this_lMieScatteringLength_BubleColumn[0], static_cast<int>(lRefractionIndex.size()))->SetSpline(true);
            MPT_holeice_temp->AddConstProperty("MIEHG_FORWARD", mMIE_spice_const[0]);
            MPT_holeice_temp->AddConstProperty("MIEHG_BACKWARD", mMIE_spice_const[1]);
            MPT_holeice_temp->AddConstProperty("MIEHG_FORWARD_RATIO", mMIE_spice_const[2]);
            mMPT_holeice_slices.push_back(MPT_holeice_temp);
            lBubleColumnMie_temp->SetMaterialPropertiesTable(mMPT_holeice_slices.at(k));
            lBubleColumnMie_slices.push_back(lBubleColumnMie_temp);
        }
    }
    //MPT_holeice_temp->DumpTable();
}

/*
 * %%%%%%%%%%%%%%%% Functions for icecube ice optical properties %%%%%%%%%%%%%%%%
 */
/**
 * This gives you temperature of ice depending on the depth.
 * Function needed for the calculation of scattering and absorption length of the ice. 
 * @param pDepth Depth in m from where we need the temperature
 * @return Temperature
 */
G4double IceCubeIce::Spice_Temperature(G4double pDepth)
{
    G4double spice_temp = 221.5 - 0.00045319 / m * pDepth + 5.822e-6 / m2 * pow(pDepth, 2.);
    return spice_temp;
}

/**
 * Calculation of the absorption length of IC-ice for a specific wavelength
 * @param pLambd Wavelength
 * @return Absorption length
 */
G4double IceCubeIce::Spice_Absorption(G4double pLambd)
{
    G4double lKappa = 1.08;
    G4double lParamA = 6954. / m;
    G4double lParamB = 6618 * nm;
    G4double lAdust = 1. / (mSpice_a400inv[mSpiceDepth_pos]) * pow(pLambd / (400. * nm), -lKappa);
    G4double lDeltaTau = Spice_Temperature(mSpice_Depth[mSpiceDepth_pos]) - Spice_Temperature(1730.);
    G4double la_inv = 1. / (lAdust + lParamA * exp(-lParamB / pLambd) * (1. + 0.01 * lDeltaTau));
    return la_inv;
}

/**
 * Calculation of the refraction index of IC-ice for a specific wavelength.
 * @param pLambd Wavelength
 * @return Refraction index
 */
G4double IceCubeIce::Spice_Refraction(G4double pLambd)
{
    // unknown depth. Parametrization by Thomas Kittler.
    G4double lLambd3 = pLambd * 1e-3;
    G4double lNphase = 1.55749 - 1.57988 / nm * lLambd3 + 3.99993 / (nm * nm) * pow(lLambd3, 2) - 4.68271 / (nm * nm * nm) * pow(lLambd3, 3) + 2.09354 / (nm * nm * nm * nm) * pow(lLambd3, 4);
    return lNphase; // using this now after discussion with Timo
}
/**
 * Calculation of the mie scattering length of IC-ice for a specific wavelength
 * @param pLambd Wavelength
 * @return Mie scattering length
 */
G4double IceCubeIce::Mie_Scattering(G4double pLambd)
{
    // depth_pos is the coordinate for the chosen depth in Depth_spice. For example to choose
    // depth=2278.2 m, we use depth_pos = 88
    G4double lAlpha = 0.90;
    G4double lAv_costheta = 0.9;
    G4double lBe_inv = 1. / (1. / (mSpice_be400inv[mSpiceDepth_pos]) * pow((pLambd / (400. * nm)), -lAlpha));
    G4double lB_inv = lBe_inv * (1. - lAv_costheta);
    return lB_inv;
}
/*
 * %%%%%%%%%%%%%%%% Functions of derived class ReflectiveSurface %%%%%%%%%%%%%%%%
 */
/**
 * Defines new reflective surface from data in json-file. 
 */
void ReflectiveSurface::ExtractInformation()
{

    pt::read_json(mFileName, mJsonTree); //read json file into mJsonTree

    mObjectName = mJsonTree.get<G4String>("jName");
    G4String lModelStr = mJsonTree.get<G4String>("jModel");
    G4String lFinishStr = mJsonTree.get<G4String>("jFinish");
    G4String lTypeStr = mJsonTree.get<G4String>("jType");

    G4OpticalSurfaceModel lModel = GetOpticalSurfaceModel(lModelStr);
    G4OpticalSurfaceFinish lFinish = GetOpticalSurfaceFinish(lFinishStr);
    G4SurfaceType lType = GetSurfaceType(lTypeStr);

    mOpticalSurface = new G4OpticalSurface(mObjectName, lModel, lFinish, lType);
    G4MaterialPropertiesTable *lMPT = new G4MaterialPropertiesTable();

    try // Only few materials have jSigmaAlpha defined
    {
        G4double lSigmaAlpha = mJsonTree.get<G4double>("jSigmaAlpha");
        mOpticalSurface->SetSigmaAlpha(lSigmaAlpha);
    }
    catch (...)
    {
    } // not very elegant, I know...

    for (pt::ptree::value_type &key : mJsonTree.get_child("jProperties"))
    {
        G4String lKey = key.second.get_value<G4String>();
        std::vector<G4double> lPhotonEnergy;
        std::vector<G4double> lValues;
        ParseToVector(lValues, mJsonTree, "jValues_" + lKey, 1., false);
        ParseToVector(lPhotonEnergy, mJsonTree, "jWavelength_" + lKey, mHC_eVnm, true);
        lMPT->AddProperty(lKey, &lPhotonEnergy[0], &lValues[0], static_cast<int>(lPhotonEnergy.size()));
    }

    mOpticalSurface->SetMaterialPropertiesTable(lMPT);
    G4String mssg = "New Optical Surface: " + mObjectName;
    info(mssg);
}

/**
 * OpticalSurfaceFinish in string to G4OpticalSurfaceFinish
 * @param  G4String
 * @return G4OpticalSurfaceFinish
 */
G4OpticalSurfaceFinish ReflectiveSurface::GetOpticalSurfaceFinish(G4String pFinish)
{
    G4OpticalSurfaceFinish lFinish;
    if (pFinish == "polished")
        lFinish = polished;
    else if (pFinish == "polishedfrontpainted")
        lFinish = polishedfrontpainted;
    else if (pFinish == "polishedbackpainted")
        lFinish = polishedbackpainted;
    else if (pFinish == "ground")
        lFinish = ground;
    else if (pFinish == "groundfrontpainted")
        lFinish = groundfrontpainted;
    else if (pFinish == "groundbackpainted")
        lFinish = groundbackpainted;
    else if (pFinish == "polishedlumirrorair")
        lFinish = polishedlumirrorair;
    else if (pFinish == "polishedlumirrorglue")
        lFinish = polishedlumirrorglue;
    else if (pFinish == "polishedair")
        lFinish = polishedair;
    else if (pFinish == "polishedteflonair")
        lFinish = polishedteflonair;
    else if (pFinish == "polishedtioair")
        lFinish = polishedtioair;
    else if (pFinish == "polishedtyvekair")
        lFinish = polishedtyvekair;
    else if (pFinish == "polishedvm2000air")
        lFinish = polishedvm2000air;
    else if (pFinish == "polishedvm2000glue")
        lFinish = polishedvm2000glue;
    else if (pFinish == "etchedlumirrorair")
        lFinish = etchedlumirrorair;
    else if (pFinish == "etchedlumirrorglue")
        lFinish = etchedlumirrorglue;
    else if (pFinish == "etchedair")
        lFinish = etchedair;
    else if (pFinish == "etchedteflonair")
        lFinish = etchedteflonair;
    else if (pFinish == "etchedtioair")
        lFinish = etchedtioair;
    else if (pFinish == "etchedtyvekair")
        lFinish = etchedtyvekair;
    else if (pFinish == "etchedvm2000air")
        lFinish = etchedvm2000air;
    else if (pFinish == "etchedvm2000glue")
        lFinish = etchedvm2000glue;
    else if (pFinish == "groundlumirrorair")
        lFinish = groundlumirrorair;
    else if (pFinish == "groundlumirrorglue")
        lFinish = groundlumirrorglue;
    else if (pFinish == "groundair")
        lFinish = groundair;
    else if (pFinish == "groundteflonair")
        lFinish = groundteflonair;
    else if (pFinish == "groundtioair")
        lFinish = groundtioair;
    else if (pFinish == "groundtyvekair")
        lFinish = groundtyvekair;
    else if (pFinish == "groundvm2000air")
        lFinish = groundvm2000air;
    else if (pFinish == "groundvm2000glue")
        lFinish = groundvm2000glue;
    return lFinish;
}
/**
 * OpticalSurfaceModel in string to G4OpticalSurfaceModel
 * @param  G4String
 * @return G4OpticalSurfaceModel
 */
G4OpticalSurfaceModel ReflectiveSurface::GetOpticalSurfaceModel(G4String pModel)
{
    G4OpticalSurfaceModel lModel;
    if (pModel == "glisur")
        lModel = glisur;
    else if (pModel == "unified")
        lModel = unified;
    else if (pModel == "LUT")
        lModel = LUT;
    return lModel;
}
/**
 * SurfaceType in string to G4SurfaceType
 * @param  G4String
 * @return G4SurfaceType
 */
G4SurfaceType ReflectiveSurface::GetSurfaceType(G4String pType)
{

    G4SurfaceType lType;
    if (pType == "dielectric_metal")
        lType = dielectric_metal;
    else if (pType == "dielectric_dielectric")
        lType = dielectric_dielectric;
    else if (pType == "dielectric_LUT")
        lType = dielectric_LUT;
    else if (pType == "firsov")
        lType = firsov;
    else if (pType == "x_ray")
        lType = x_ray;
    return lType;
}

