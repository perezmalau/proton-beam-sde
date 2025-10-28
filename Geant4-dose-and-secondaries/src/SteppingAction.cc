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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "Analysis.hh"

SteppingAction::SteppingAction(TrackingAction* trackingAction)
  : G4UserSteppingAction(),
    fTrackingAction(trackingAction) {}


SteppingAction::~SteppingAction() {}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // G4StepPoint* prePoint  = step->GetPreStepPoint();
    // G4StepPoint* postPoint = step->GetPostStepPoint();
    //
    // G4VPhysicalVolume* preVolume  = prePoint->GetPhysicalVolume();
    // G4VPhysicalVolume* postVolume = postPoint->GetPhysicalVolume();
    // auto parentID = step->GetTrack()->GetParentID();
    //
    // auto init_energy = fTrackingAction->GetInitEnergy();
    // auto posX = fTrackingAction->GetInitPosX();
    // auto posY = fTrackingAction->GetInitPosY();
    // auto posZ = fTrackingAction->GetInitPosZ();
    // auto eloss = init_energy - postPoint->GetKineticEnergy();
    // auto eloss_step = prePoint->GetKineticEnergy() - postPoint->GetKineticEnergy();
    //
    // G4String process = postPoint->GetProcessDefinedStep()->GetProcessName();
    // G4String particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    // G4int evtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //
    // // Check up on hadron elastic processes
    // if (process=="hadElastic" && particle=="proton" && parentID==0) {
    //     G4ThreeVector preMomentum = prePoint->GetMomentumDirection();
    //     G4ThreeVector postMomentum = postPoint->GetMomentumDirection();
    //     // Both are unit vectors, so the dot product gives cos(theta)
    //     double cosTheta = preMomentum.dot(postMomentum);
    //
    //     // Clamp value to valid range for acos to avoid numerical errors
    //     if (cosTheta > 1.0) cosTheta = 1.0;
    //     if (cosTheta < -1.0) cosTheta = -1.0;
    //
    //     auto analysisManager = G4AnalysisManager::Instance();
    //     analysisManager->FillNtupleIColumn(1, 0, evtID);
    //     analysisManager->FillNtupleDColumn(1, 1, prePoint->GetKineticEnergy());
    //     analysisManager->FillNtupleDColumn(1, 2, eloss_step);
    //     analysisManager->FillNtupleDColumn(1, 3, cosTheta);
    //     analysisManager->AddNtupleRow(1);
    // }
    //
    // // Check up on Proton Inelastic processes
    // G4ThreeVector vertex = step->GetPostStepPoint()->GetPosition();
    // if (process=="protonInelastic" && parentID==0) {
    //     auto protonEnergy = prePoint->GetKineticEnergy();
    //     const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
    //     for (auto sec : *secondaries) {
    //         G4String name = sec->GetDefinition()->GetParticleName();
    //         G4double energy = sec->GetKineticEnergy();
    //         auto analysisManager = G4AnalysisManager::Instance();
    //         analysisManager->FillNtupleIColumn(2, 0, evtID);
    //         analysisManager->FillNtupleSColumn(2, 1, name);
    //         analysisManager->FillNtupleDColumn(2, 2, energy);
    //         analysisManager->FillNtupleDColumn(2, 3, vertex.x());
    //         analysisManager->FillNtupleDColumn(2, 4, vertex.y());
    //         analysisManager->FillNtupleDColumn(2, 5, vertex.z());
    //         analysisManager->FillNtupleDColumn(2, 6, protonEnergy);
    //         analysisManager->AddNtupleRow(2);
    //     }
    // }
    //
    // // For those secondary particles that managed to escape the phantom, record creation positions and energy left in phantom
    // if (preVolume->GetName() == "voxel" && postVolume->GetName() == "World" && parentID > 0) {
    //     auto analysisManager = G4AnalysisManager::Instance();
    //     analysisManager->FillNtupleIColumn(0, 0, evtID);
    //     analysisManager->FillNtupleIColumn(0, 1, parentID);
    //     analysisManager->FillNtupleSColumn(0, 2, particle);
    //     analysisManager->FillNtupleDColumn(0, 3, init_energy);
    //     analysisManager->FillNtupleDColumn(0, 4, eloss);
    //     analysisManager->FillNtupleDColumn(0, 5, posX);
    //     analysisManager->FillNtupleDColumn(0, 6, posY);
    //     analysisManager->FillNtupleDColumn(0, 7, posZ);
    //     analysisManager->AddNtupleRow(0);
    // }
}
