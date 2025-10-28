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
/// \file src/RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4Timer.hh"
#include "Analysis.hh"
#include "G4RunManager.hh"
#include "MyRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction(),
	fDoseOutputFileName("edep_map.csv"),
	fRootOutputFileName("ExtraInfoRun.root")
{
	fMessenger = new G4GenericMessenger(this, "/output/", "Output control");
	fMessenger->DeclareProperty("doseFileName", fDoseOutputFileName, "Set dose output csv file name");
	fMessenger->DeclareProperty("rootFileName", fRootOutputFileName, "Set root output file name");

    timer = new G4Timer();
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetFileName(fRootOutputFileName);
	G4cout << "Using " << analysisManager->GetType() << G4endl;
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetNtupleMerging(true); //only with root

	// Creating ntuple
	//
	// analysisManager->CreateNtuple("SecondaryParticles", "SecondaryParticles");
	// analysisManager->CreateNtupleIColumn("EventID");
	// analysisManager->CreateNtupleIColumn("ParentID");
	// analysisManager->CreateNtupleSColumn("Particle");
	// analysisManager->CreateNtupleDColumn("InitEnergy");
	// analysisManager->CreateNtupleDColumn("deltaE");
	// analysisManager->CreateNtupleDColumn("InitPosX");
	// analysisManager->CreateNtupleDColumn("InitPosY");
	// analysisManager->CreateNtupleDColumn("InitPosZ");
	// analysisManager->FinishNtuple();
	//
	// analysisManager->CreateNtuple("hadElasticData", "hadElasticData");
	// analysisManager->CreateNtupleIColumn("EventID");
	// analysisManager->CreateNtupleDColumn("preStepEnergy");
	// analysisManager->CreateNtupleDColumn("eneTransfer");
	// analysisManager->CreateNtupleDColumn("CosTheta");
	// analysisManager->FinishNtuple();
	//
	// analysisManager->CreateNtuple("protonInelasticData", "protonInelasticData");
	// analysisManager->CreateNtupleIColumn("EventID");
	// analysisManager->CreateNtupleSColumn("Particle");
	// analysisManager->CreateNtupleDColumn("KineticEnergy");
	// analysisManager->CreateNtupleDColumn("posX");
	// analysisManager->CreateNtupleDColumn("posY");
	// analysisManager->CreateNtupleDColumn("posZ");
	// analysisManager->CreateNtupleDColumn("InitProtonEnergy");
	// analysisManager->FinishNtuple();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
	delete G4AnalysisManager::Instance();
	delete fMessenger;
}


G4Run* RunAction::GenerateRun()
{
	return new MyRun();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	if (isMaster)
        	timer->Start();
	// Get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->OpenFile(fRootOutputFileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  	if (isMaster) {
		timer->Stop();
	  	G4cout << "Simulation time: " << timer->GetRealElapsed()<< " seconds." << G4endl;
	}

	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

	const MyRun* myRun = static_cast<const MyRun*>(run);
	const auto& edepMap = myRun->fTotalEdep;

	std::ofstream outFile(fDoseOutputFileName);
	outFile << "voxelX,voxelY,voxelZ,eDep\n";

	for (const auto& [key, edep] : edepMap) {
		int i, j, k;
		std::tie(i, j, k) = key;
		outFile << i << "," << j << "," << k << "," << edep << "\n";
	}

	outFile.close();
	G4cout << "Energy deposition map saved to " << fDoseOutputFileName << "\n";
}
