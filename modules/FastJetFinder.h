/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FastJetFinder_h
#define FastJetFinder_h

/** \class FastJetFinder
 *
 *  Finds jets using FastJet library.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <memory>
#include <vector>

#include "fastjet/JetDefinition.hh"

#include "TObjArray.h"

class TIterator;

#include "fastjet/PseudoJet.hh"

namespace fastjet
{
class AreaDefinition;
class JetMedianBackgroundEstimator;
namespace contrib
{
class NjettinessPlugin;
class ValenciaPlugin;
class AxesDefinition;
class MeasureDefinition;
} // namespace contrib
} // namespace fastjet

class FastJetFinder: public DelphesModule
{
public:
  FastJetFinder();
  ~FastJetFinder();

  void Init();
  void Process();
  void Finish();

private:
  std::unique_ptr<fastjet::JetDefinition::Plugin> fPlugin; //!
  std::unique_ptr<fastjet::JetDefinition::Recombiner> fRecomb; //!

  std::unique_ptr<fastjet::contrib::AxesDefinition> fAxesDef; //!
  std::unique_ptr<fastjet::contrib::MeasureDefinition> fMeasureDef; //!

  std::unique_ptr<fastjet::contrib::NjettinessPlugin> fNjettinessPlugin; //!
  std::unique_ptr<fastjet::contrib::ValenciaPlugin> fValenciaPlugin; //!
  std::unique_ptr<fastjet::JetDefinition> fDefinition; //!

  Int_t fJetAlgorithm;
  Double_t fParameterR;
  Double_t fParameterP;

  Double_t fJetPTMin;
  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Int_t fMaxIterations;
  Int_t fMaxPairSize;
  Int_t fIratch;
  Int_t fAdjacencyCut;
  Double_t fOverlapThreshold;

  //-- Exclusive clustering for e+e- collisions --

  Int_t fNJets;
  Double_t fDCut;
  Bool_t fExclusiveClustering;

  //-- Valencia Linear Collider algorithm
  Double_t fGamma;

  //-- N (sub)jettiness parameters --

  Bool_t fComputeNsubjettiness;
  Double_t fBeta;
  Int_t fAxisMode;
  Double_t fRcutOff;
  Int_t fN;

  //-- Trimming parameters --

  Bool_t fComputeTrimming;
  Double_t fRTrim;
  Double_t fPtFracTrim;

  //-- Pruning parameters --

  Bool_t fComputePruning;
  Double_t fZcutPrun;
  Double_t fRcutPrun;
  Double_t fRPrun;

  //-- SoftDrop parameters --

  Bool_t fComputeSoftDrop;
  Double_t fBetaSoftDrop;
  Double_t fSymmetryCutSoftDrop;
  Double_t fR0SoftDrop;

  // --- FastJet Area method --------

  std::unique_ptr<fastjet::AreaDefinition> fAreaDefinition; //!
  Int_t fAreaAlgorithm;
  Bool_t fComputeRho;

  // -- ghost based areas --
  Double_t fGhostEtaMax;
  Int_t fRepeat;
  Double_t fGhostArea;
  Double_t fGridScatter;
  Double_t fPtScatter;
  Double_t fMeanGhostPt;

  // -- voronoi areas --
  Double_t fEffectiveRfact;

  struct TEntryStruct
  {
    std::unique_ptr<fastjet::JetMedianBackgroundEstimator> estimator;
    Double_t etaMin, etaMax;
  };

  std::vector<TEntryStruct> fEstimators; //!

  std::unique_ptr<TIterator> fItInputArray; //!

  const TObjArray *fInputArray = nullptr; //!

  TObjArray *fOutputArray = nullptr; //!
  TObjArray *fRhoOutputArray = nullptr; //!
  TObjArray *fConstituentsOutputArray = nullptr; //!

  // ---------------------------------------------------------------------------
  // Dark-hadron visible-jet matching  (algorithm ID kDarkHadronVisibleMatch)
  //
  // Inputs (configured in the card):
  //   InputArray         : visible final-state SM particles (from PdgCodeFilter)
  //   DarkHadronJetArray : pre-clustered dark-hadron jets   (from FastJetFinder)
  //   ParticleInputArray : full GenParticle collection      (from Delphes reader)
  //
  // For each dark-hadron jet the algorithm:
  //   1. Collects the GenParticle indices of the jet's DH constituents.
  //   2. BFS-traces every visible SM particle's ancestry (via M1/M2) to find
  //      the owning DH jet.  Results are memoised per GenParticle index so
  //      each node is visited at most once per event.
  //   3. Builds one output PseudoJet per DH jet via fastjet::join(), so
  //      jet.constituents() is populated for all post-clustering tools.
  //      Output ordering mirrors fDarkHadronJetArray (NOT re-sorted by pT).
  // ---------------------------------------------------------------------------

  static const Int_t kDarkHadronVisibleMatch = 20;

  /// Pre-clustered dark-hadron jets (already anti-kT clustered).
  TObjArray *fDarkHadronJetArray = nullptr;

  /// Full GenParticle collection for ancestry tracing.
  TObjArray *fAllParticleArray = nullptr;

  /// Per-event working array: visible Candidates matched to jets,
  /// stored in user_index order.  Owned by the factory (nodelete).
  TObjArray fDHMatchedVisibleArray;

  /// These are heap-allocated per-event and deleted at the start of the
  /// following event, after ROOT has written the previous event's jets.
  std::vector<fastjet::ClusterSequence *> fDHClusterSequences;

  /// Core implementation: fills outputJets and matchedVisArray.
  void BuildDarkHadronMatchedJets(
    std::vector<fastjet::PseudoJet> &outputJets,
    TObjArray                       &matchedVisArray);

  ClassDef(FastJetFinder, 1)
};

#endif
