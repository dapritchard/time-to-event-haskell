{-# LANGUAGE DuplicateRecordFields #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}
-- |

module CoxPH.CoxPHUpdateNewtonRaphson where

import CoxPH.Data
import Numeric.LinearAlgebra ( Matrix, add, cols, diag, outer, scale )
import Numeric.LinearAlgebra.Data qualified as LD
import Data.Vector.Storable qualified as VS
import Data.Vector qualified as V

data IterationInfo = IterationInfo
  { subjectIndex :: Int
  , time :: Double
  , stratum :: Int
  , nEvents :: Int
  }
  deriving Show

data TTEData = TTEData
  { time :: VS.Vector Double
  , eventStatus :: V.Vector Delta
  , xDesignMatrix :: Matrix Double
  , stratum :: VS.Vector Int
  , weights :: VS.Vector Double
  , xProdBeta :: VS.Vector Double
  , tiesMethod :: CoxPHMethod
  }
  deriving Show

data NRTerms = NRTerms
  { sumWeights :: Double
  , sumWeightedRisk :: Double
  , logLikelihood :: Double
  , score :: VS.Vector Double
  , xBarUnscaled :: VS.Vector Double
  , informationTerm1 :: Matrix Double
  }
  deriving Show

data NRResults = NRResults
  { sumLogLikelihood :: Double
  , score :: VS.Vector Double
  , informationMatrix :: Matrix Double
  }
  deriving Show

coxPHUpdateNewtonRaphson :: TTEData -> NRResults
coxPHUpdateNewtonRaphson tteData =
  let p = cols tteData.xDesignMatrix
      iterationInfo = IterationInfo
        { subjectIndex = VS.length tteData.time - 1
        , time = 0
        , stratum = 0
        , nEvents = 0
        }
      nrResults = NRResults
        { sumLogLikelihood = 0
        , score = VS.replicate p 0
        , informationMatrix = createEmptyMatrix p
        }
      (_, newNRResults) = calcStrata tteData iterationInfo nrResults
  in  newNRResults

calcStrata
  :: TTEData
  -> IterationInfo
  -> NRResults
  -> (IterationInfo, NRResults)
calcStrata tteData iterationInfo nrResults
  -- Case: we've seen all of the subjects. This is the base case
  | iterationInfo.subjectIndex < 0 = -- TODO: create helper functions for checks
    (iterationInfo, nrResults)
  -- Case: compute the remaining strata. We calculate all of the
  -- relevant terms for the current stratum  and recursively call `calcStrata` to
  -- conditionally compute the next stratum
  | otherwise =
    let p = VS.length nrResults.score
        initialIterationInfo = resetIterationInfo iterationInfo
        initialOverallData = createInitialData p
        initialInformation = createEmptyMatrix p
        (newIterationInfo, newNRTerms, informationMatrix) =
          calcTimeBlocks tteData
                         initialIterationInfo
                         initialOverallData
                         initialInformation
        aggregatedNRResults = aggregateNRResults nrResults
                                                 newNRTerms
                                                 informationMatrix
    in calcStrata tteData newIterationInfo aggregatedNRResults
  where
    resetIterationInfo :: IterationInfo -> IterationInfo
    resetIterationInfo iterationInfo =
      IterationInfo
        { subjectIndex = iterationInfo.subjectIndex
        , time = iterationInfo.time
        , stratum = tteData.stratum VS.! iterationInfo.subjectIndex
        , nEvents = iterationInfo.nEvents
        }
    aggregateNRResults
      :: NRResults
      -> NRTerms
      -> Matrix Double
      -> NRResults
    aggregateNRResults nrResults overallData informationMatrix =
      NRResults
        { sumLogLikelihood = nrResults.sumLogLikelihood
                             + overallData.logLikelihood
        , score = add nrResults.score overallData.score
        , informationMatrix = add nrResults.informationMatrix informationMatrix
        }

calcTimeBlocks
  :: TTEData
  -> IterationInfo
  -> NRTerms
  -> Matrix Double
  -> (IterationInfo, NRTerms, Matrix Double)
calcTimeBlocks tteData iterationInfo overallData informationMatrix
  -- Case: we've either seen all of the subjects or we've found a subject with a
  -- different stratum. This is the base case
  | (iterationInfo.subjectIndex < 0) -- TODO: create helper functions for checks
      || checkDifferentStrata tteData iterationInfo =
      (iterationInfo, overallData, informationMatrix)
  -- Case: the current subject is part of the current stratum (note that it
  -- could be the first subject in the stratum). We calculate and accumulate all
  -- of the relevant terms for the subjects within a strata who were censored or
  -- had an event at the current time, and recursively call
  -- `calcTimeBlocksSubjects` to conditionally compute the next time block
  | otherwise =
      let p = VS.length overallData.score
          initialTiedData = createInitialData p
          initialIterationInfo = resetIterationInfo iterationInfo
          (timeBlockIterationInfo, timeBlockNRTerms, timeBlockTiedData) =
            calcTimeBlocksSubjects tteData
                                   initialIterationInfo
                                   overallData
                                   initialTiedData
          aggregateData = case tteData.tiesMethod of
            Breslow -> computeBreslow
            Efron   -> computeEfron
          (newNRTerms, newInformationMatrix) =
            aggregateData timeBlockIterationInfo
                          timeBlockNRTerms
                          timeBlockTiedData
                          informationMatrix
      in  calcTimeBlocks tteData
                         timeBlockIterationInfo
                         newNRTerms
                         newInformationMatrix
  where
    resetIterationInfo :: IterationInfo -> IterationInfo
    resetIterationInfo iterationInfo =
      IterationInfo
        { subjectIndex = iterationInfo.subjectIndex
        , time = tteData.time VS.! iterationInfo.subjectIndex
        , stratum = iterationInfo.stratum
        , nEvents = 0
        }

computeBreslow
  :: IterationInfo
  -> NRTerms
  -> NRTerms
  -> Matrix Double
  -> (NRTerms, Matrix Double)
computeBreslow iterationInfo overallData tiedData informationMatrix
  | iterationInfo.nEvents == 0 =
    let newLogLikelihood = overallData.logLikelihood + tiedData.logLikelihood -- TODO: is there any need to update @overallData@ at all?
        newNRTerms = (overallData :: NRTerms) { logLikelihood = newLogLikelihood }
    in  (newNRTerms, informationMatrix)
  | otherwise =
    let
        newSumWeightedRisk = overallData.sumWeightedRisk
                             + tiedData.sumWeightedRisk
        newXBarUnscaled = add overallData.xBarUnscaled
                              tiedData.xBarUnscaled
        newXBar = scale (1 / newSumWeightedRisk) newXBarUnscaled
        newInformationTerm1 = add overallData.informationTerm1
                                  tiedData.informationTerm1
        newNRTerms = NRTerms
          { sumWeights = 0
          , sumWeightedRisk = newSumWeightedRisk
          , logLikelihood = overallData.logLikelihood
                            + tiedData.logLikelihood
                            - (tiedData.sumWeights * log newSumWeightedRisk)
          , score = add overallData.score
                        (add tiedData.score
                             (scale (- tiedData.sumWeights) newXBar))
          , xBarUnscaled = newXBarUnscaled
          , informationTerm1 = newInformationTerm1
          }
        blockInformation = scale (tiedData.sumWeights / newSumWeightedRisk)
                                 (add newInformationTerm1
                                      (scale (- newSumWeightedRisk)
                                             (outer newXBar newXBar)))
        newInformationMatrix = add informationMatrix blockInformation
    in  (newNRTerms, newInformationMatrix)

computeEfron
  :: IterationInfo
  -> NRTerms
  -> NRTerms
  -> Matrix Double
  -> (NRTerms, Matrix Double)
computeEfron iterationInfo overallData tiedData informationMatrix
  | iterationInfo.nEvents <= 1 =
    computeBreslow iterationInfo overallData tiedData informationMatrix
  | otherwise =
    undefined -- FIXME

calcTimeBlocksSubjects
  :: TTEData
  -> IterationInfo
  -> NRTerms
  -> NRTerms
  -> (IterationInfo, NRTerms, NRTerms)
calcTimeBlocksSubjects tteData iterationInfo overallData tiedData
  -- Case: we've either seen all of the subjects or we've found a subject with a
  -- different censoring or event time. This is the base case
  | (iterationInfo.subjectIndex < 0) -- TODO: create helper functions for checks
      || checkDifferentStrata tteData iterationInfo
      || (tteData.time VS.! iterationInfo.subjectIndex) /= iterationInfo.time =
      (iterationInfo, overallData, tiedData)
  -- Case: the current subject is part of the current censoring or event tied
  -- time block (note that it could be the first subject in the block). We
  -- calculate and accumulate all of the relevant terms for the subject and
  -- recursively call `calcTimeBlocksSubjects` to conditionally compute the next
  -- subject in the time block
  | otherwise =
      let subjectWeight = tteData.weights VS.! iterationInfo.subjectIndex
          subjectXProdBeta = tteData.xProdBeta VS.! iterationInfo.subjectIndex
          subjectWeightedRisk = subjectWeight * exp subjectXProdBeta
          subjectXRow = tteData.xDesignMatrix LD.! iterationInfo.subjectIndex
          subjectWeightedXRow = scale subjectWeightedRisk subjectXRow
          newInformationTerm1 = scale subjectWeightedRisk
                                      (outer subjectXRow subjectXRow)
      in  case tteData.eventStatus V.! iterationInfo.subjectIndex of
            Censored ->
              let newOverallData = overallData
                    { sumWeightedRisk = overallData.sumWeightedRisk
                                        + subjectWeightedRisk
                    , xBarUnscaled = add overallData.xBarUnscaled
                                         subjectWeightedXRow
                    , informationTerm1 = add overallData.informationTerm1
                                             newInformationTerm1
                    }
                  newIterationInfo = iterationInfo
                    { subjectIndex = iterationInfo.subjectIndex - 1
                    }
              in  calcTimeBlocksSubjects
                    tteData
                    newIterationInfo
                    newOverallData
                    tiedData
            ObservedEvent ->
              let newTiedData = tiedData
                    { sumWeights = tiedData.sumWeights + subjectWeight
                    , sumWeightedRisk = tiedData.sumWeights
                                        + subjectWeightedRisk
                    , logLikelihood = tiedData.logLikelihood
                                      + (subjectWeight * subjectXProdBeta) -- TODO: split this into logLikelihoodTerm1 and logLikelihoodTerm2?
                    , score = add tiedData.score (scale subjectWeight subjectXRow)
                    , xBarUnscaled = add tiedData.xBarUnscaled
                                         subjectWeightedXRow
                    , informationTerm1 = newInformationTerm1
                    }
                  newIterationInfo = iterationInfo
                    { subjectIndex = iterationInfo.subjectIndex - 1
                    , nEvents = iterationInfo.nEvents + 1
                    }
              in  calcTimeBlocksSubjects
                    tteData
                    newIterationInfo
                    overallData
                    newTiedData

createInitialData :: Int -> NRTerms
createInitialData p
  = NRTerms
      { sumWeights = 0
      , sumWeightedRisk = 0
      , logLikelihood = 0
      , score = VS.replicate p 0
      , xBarUnscaled = VS.replicate p 0
      , informationTerm1 = createEmptyMatrix p
      }

createEmptyMatrix :: Int -> Matrix Double
createEmptyMatrix p = diag (VS.replicate p 0)

checkDifferentStrata :: TTEData -> IterationInfo -> Bool
checkDifferentStrata tteData iterationInfo =
  let subjectStratum = tteData.stratum VS.! iterationInfo.subjectIndex
  in  subjectStratum /= iterationInfo.stratum

-- v1 :: Vector Double
-- v1 = L.fromList [3, 6]

-- v2 :: Vector Double
-- v2 = L.fromList [3, 6]

-- m :: Matrix Double
-- m = fromColumns [v1, v2]

-- testResult :: Vector Double
-- testResult =  testMatrix #> testVector
