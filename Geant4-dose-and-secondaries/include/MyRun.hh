#include <map>
#include <tuple>
#include "G4Run.hh"
#pragma once

class MyRun : public G4Run {
public:
    std::map<std::tuple<int, int, int>, double> fTotalEdep;

    void AddEdep(int x, int y, int z, double edep) {
        fTotalEdep[std::make_tuple(x, y, z)] += edep;
    }

    void Merge(const G4Run* localRun) override {
        const auto* local = dynamic_cast<const MyRun*>(localRun);
        for (const auto& [key, val] : local->fTotalEdep) {
            fTotalEdep[key] += val;
        }
        G4Run::Merge(localRun);
    }
};