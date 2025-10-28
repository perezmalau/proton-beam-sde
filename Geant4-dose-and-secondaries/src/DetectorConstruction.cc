//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorSD.hh"
#include "G4UserLimits.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"

#include "G4ScoringManager.hh"
#include "G4ScoringBox.hh"
#include "G4PSEnergyDeposit.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),
	fCheckOverlaps(true),
    fPhantomSize(20.*cm),
    fAirThickness(0.*cm),
    fBoneThickness(0.*cm)
{
    fMessenger = new G4GenericMessenger(this, "/phantom/", "Phantom geometry control");
    fMessenger->DeclarePropertyWithUnit("phantomSize", "cm",
        fPhantomSize, "Set initial size of cubic water phantom.");
    fMessenger->DeclarePropertyWithUnit("airThickness", "cm",
        fAirThickness, "Set initial air gap thickness to be added.");
    fMessenger->DeclarePropertyWithUnit("boneThickness", "cm",
    fBoneThickness, "Set bone slab thickness to be added.");
    DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    //clean out old geometry before updating
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    // Define volumes
    return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    // standard material defined using NIST Manager
    G4NistManager *nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_WATER");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_B-100_BONE");
	G4String name, symbol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    auto phantomSD = new DetectorSD("PhantomSD");
    G4SDManager::GetSDMpointer()->AddNewDetector(phantomSD);
    SetSensitiveDetector("water", phantomSD);
    if (fAirThickness > 0.){
        SetSensitiveDetector("air", phantomSD);
    }
    if (fBoneThickness > 0.){
        SetSensitiveDetector("bone", phantomSD);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCubicPhantomSize(G4double val) {
    fPhantomSize = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBoneThickness(G4double val) {
    if (val + fAirThickness > fPhantomSize) {
        G4ExceptionDescription msg;
        msg << "Invalid geometry: bone thickness (" << val/cm
            << " cm) + air thickness (" << fAirThickness/cm
        << " cm) exceeds phantom size (" << fPhantomSize/cm
        << " cm)";
        G4Exception("DetectorConstruction::SetBoneThickness",
                    "Geom001", JustWarning, msg);
        return; // reject the new value
    }
    fBoneThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAirThickness(G4double val) {
    if (val + fBoneThickness > fPhantomSize) {
        G4ExceptionDescription msg;
        msg << "Invalid geometry: air thickness (" << val/cm
            << " cm) + bone thickness (" << fBoneThickness/cm
        << " cm) exceeds phantom size (" << fPhantomSize/cm
    << " cm)";
        G4Exception("DetectorConstruction::SetAirThickness",
                    "Geom002", JustWarning, msg);
        return; // reject the new value
    }
    fAirThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{

    // Geometry parameters
    G4double worldSizeXY = fPhantomSize*5;
    G4double worldSizeZ = fPhantomSize*5;

    // Get materials - change them if needed
    G4Material *defaultMaterial = G4Material::GetMaterial("G4_AIR");
	G4Material *tissueMaterial = G4Material::GetMaterial("G4_WATER");
    G4Material *boneMaterial = G4Material::GetMaterial("G4_B-100_BONE");

    //
    // WORLD
    //
    G4VSolid *worldSolid
        = new G4Box("World", worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

    auto *worldLogical = new G4LogicalVolume(worldSolid, defaultMaterial, "World");       // its name

    G4VPhysicalVolume *worldPV = new G4PVPlacement(nullptr,        // no rotation
            G4ThreeVector(),  // at (0,0,0)
            worldLogical,          // its logical volume
            "World",          // its name
            nullptr,                // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps
    //
    // PHANTOM
    //
    auto* phantSolid   = new G4Box("phantom", fPhantomSize/2., fPhantomSize/2., fPhantomSize/2.);
    auto* phantLogical = new G4LogicalVolume(phantSolid, defaultMaterial, "phantom");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), phantLogical, "phantom",
        worldLogical, false, 0);

    // Set step limit in phantom
    auto stepLimit = new G4UserLimits(0.5*mm);

    phantLogical->SetUserLimits(stepLimit);
    // Air slab
    if (fAirThickness > 0.) {
        auto* airSolid   = new G4Box("air", fAirThickness/2., fPhantomSize/2., fPhantomSize/2.);
        auto* airLogical = new G4LogicalVolume(airSolid, defaultMaterial, "air");
        new G4PVPlacement(nullptr, G4ThreeVector(-fPhantomSize/2 + fAirThickness/2, 0, 0),
            airLogical, "air",
            phantLogical, false, 0);

        auto *airVisAtt = new G4VisAttributes(G4Colour(1, 1, 1));
        airVisAtt->SetVisibility(true);
        airLogical->SetVisAttributes(airVisAtt);
        airLogical->SetUserLimits(stepLimit);

        // create region and attach logical volume
        G4Region* airRegion = new G4Region("AirRegion");
        airRegion->AddRootLogicalVolume(airLogical);

        // create production cuts and set a large range for electrons
        G4ProductionCuts* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(1.*mm, "e-");
        // assign cuts to region
        airRegion->SetProductionCuts(cuts);
    }
    // Bone slab
    if (fBoneThickness > 0.) {
        auto* boneSolid   = new G4Box("bone", fBoneThickness/2., fPhantomSize/2., fPhantomSize/2.);
        auto* boneLogical = new G4LogicalVolume(boneSolid, boneMaterial, "bone");
        new G4PVPlacement(nullptr, G4ThreeVector(-fPhantomSize/2 + fAirThickness + fBoneThickness/2, 0, 0),
            boneLogical, "bone",
            phantLogical, false, 0);

        boneLogical->SetUserLimits(stepLimit);

        auto *boneVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
        boneVisAtt->SetVisibility(true);
        boneLogical->SetVisAttributes(boneVisAtt);
    }
    // Water slab - cover the remaining bit of the phantom
    G4double tissueThickness = fPhantomSize-fAirThickness-fBoneThickness;
    auto* waterSolid   = new G4Box("water", tissueThickness/2., fPhantomSize/2., fPhantomSize/2.);
    auto* waterLogical = new G4LogicalVolume(waterSolid, tissueMaterial, "water");
    new G4PVPlacement(nullptr, G4ThreeVector(fPhantomSize/2 - tissueThickness/2, 0, 0),
        waterLogical, "water",
        phantLogical, false, 0);
    waterLogical->SetUserLimits(stepLimit);

    auto *waterVisAtt = new G4VisAttributes(G4Colour(0.2, 0.6, 0.9));
    waterVisAtt->SetVisibility(true);
    waterVisAtt->SetDaughtersInvisible(false);
    waterLogical->SetVisAttributes(waterVisAtt);


    //
    // Visualization attributes
    //
    auto *worldVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    worldVisAtt->SetVisibility(false);
    worldLogical->SetVisAttributes(worldVisAtt);


    //
    // Always return the physical World
    //
    return worldPV;
}
