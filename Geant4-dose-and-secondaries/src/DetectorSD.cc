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
/// \file src/DetectorSD.cc

#include "DetectorSD.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "MyRun.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorSD::DetectorSD(const G4String &name)
    : G4VSensitiveDetector(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorSD::~DetectorSD()
= default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DetectorSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{
    G4double edep = step->GetTotalEnergyDeposit()/MeV;
    if (edep == 0.) return false;

    G4double pixelSize = 0.1 *cm;
    G4ThreeVector voxelGridOrigin = G4ThreeVector(-10*cm, -10*cm, -10*cm);

    G4ThreeVector pos = 0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());
    G4ThreeVector voxelIndex = (pos - voxelGridOrigin) / pixelSize;

    G4int ix = std::floor(voxelIndex.x());
    G4int iy = std::floor(voxelIndex.y());
    G4int iz = std::floor(voxelIndex.z());

    if (ix < 0 || ix >= 200 ||
    iy < 0 || iy >= 200 ||
    iz < 0 || iz >= 200) {
        return false;
    }

    //G4String particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    // G4int trackID = step->GetTrack()->GetTrackID();
    MyRun* run = static_cast<MyRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    //if (particle=="proton" and trackID==1) {
    run->AddEdep(ix, iy, iz, edep);
    //}

    return true;
}
