{-# LANGUAGE DuplicateRecordFields #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}
-- |

module CoxPH.CoxPHUpdateNewtonRaphson where

import CoxPH.Data
-- import Numeric.LinearAlgebra.Data -- (Matrix, Vector, fromColumns)
import Numeric.LinearAlgebra qualified as L -- (Matrix, Vector, fromColumns)
import Numeric.LinearAlgebra ( (!), (#>), Matrix, Vector, add, diag, outer, scale, fromColumns )
import Data.Vector.Storable qualified as VS
import Data.Vector qualified as V

updateNewtonRaphson :: [StrataData] -> VS.Vector Double -> Vector Double -> ()
updateNewtonRaphson strataDatas weights beta =
  -- let z = add ()
  ()
-- updateTime

-- data TTEData = TTEData {

--   }

data IterationInfo = IterationInfo
  { subjectIndex :: Int
  , time :: Double
  , stratum :: Int
  , nEvents :: Int
  }

data StrataData = StrataData
  { time :: VS.Vector Double
  , eventStatus :: V.Vector Delta
  , xDesignMatrix :: Matrix Double
  , xOffset :: Vector Double
  , stratum :: Vector Int
  , weights :: Vector Double
  , xProdBeta :: Vector Double
  , tiesMethod :: CoxPHMethod
  }

data OverallData = OverallData
  { sumWeights :: Double
  , sumWeightedRisk :: Double
  , logLikelihood :: Double
  , score :: Vector Double
  , xBarUnscaled :: Vector Double
  , informationTerm1 :: Matrix Double
  }

-- data TiedData = TiedData
--   { sumWeights :: Double
--   , sumWeightedRisk :: Double
--   , logLikelihood :: Double
--   , score :: Vector Double
--   , xBarUnscaled :: Vector Double
--   , informationTerm1 :: Matrix Double
--   }

data NRUpdateResults = NRUpdateResults
  { sumLogLikelihood :: Double
  , score :: Vector Double
  , informationMatrix :: Matrix Double
  }

calcStrata
  :: StrataData
  -> IterationInfo
  -> NRUpdateResults
  -> (IterationInfo, NRUpdateResults)
calcStrata strataData iterationInfo nRUpdateResults
  -- Case: we've seen all of the subjects. This is the base case
  | iterationInfo.subjectIndex < 0 = -- TODO: create helper functions for checks
    (iterationInfo, nRUpdateResults)
  -- Case: compute the remaining strata. We calculate all of the
  -- relevant terms for the current stratum  and recursively call `calcStrata` to
  -- conditionally compute the next stratum
  | otherwise =
    let p = VS.length nRUpdateResults.score
        initialData = createInitialData p
        initialInformation = createEmptyMatrix p
        (newIterationInfo, newOverallData, informationMatrix) =
          calcTimeBlocks strataData iterationInfo initialData initialInformation
        aggregatedNRUpdateResults = aggregateNRUpdateResults nRUpdateResults
                                                             newOverallData
                                                             informationMatrix
    in calcStrata strataData newIterationInfo aggregatedNRUpdateResults
  where
    aggregateNRUpdateResults
      :: NRUpdateResults
      -> OverallData
      -> Matrix Double
      -> NRUpdateResults
    aggregateNRUpdateResults nRUpdateResults overallData informationMatrix =
      NRUpdateResults
        { sumLogLikelihood = nRUpdateResults.sumLogLikelihood
                             + overallData.logLikelihood
        , score = add nRUpdateResults.score overallData.score
        , informationMatrix = add nRUpdateResults.informationMatrix
                                  informationMatrix
        }

calcTimeBlocks
  :: StrataData
  -> IterationInfo
  -> OverallData
  -> Matrix Double
  -> (IterationInfo, OverallData, Matrix Double)
calcTimeBlocks strataData iterationInfo overallData informationMatrix
  -- Case: we've either seen all of the subjects or we've found a subject with a
  -- different stratum. This is the base case
  | (iterationInfo.subjectIndex < 0) -- TODO: create helper functions for checks
      || checkDifferentStrata strataData iterationInfo =
      (iterationInfo, overallData, informationMatrix)
  -- Case: the current subject is part of the current stratum (note that it
  -- could be the first subject in the stratum). We calculate and accumulate all
  -- of the relevant terms for the subjects within a strata who were censored or
  -- had an event at the current time, and recursively call
  -- `calcTimeBlocksSubjects` to conditionally compute the next time block
  | otherwise =
      let p = VS.length overallData.score
          initialTiedData = createInitialData p
          (timeBlockIterationInfo, timeBlockOverallData, timeBlockTiedData) =
            calcTimeBlocksSubjects strataData
                                   iterationInfo
                                   overallData
                                   initialTiedData
          aggregateData = case strataData.tiesMethod of
            Breslow -> computeBreslow
            Efron   -> computeEfron
          (newOverallData, newInformationMatrix) =
            aggregateData timeBlockIterationInfo
                          timeBlockOverallData
                          timeBlockTiedData
                          informationMatrix
      in  calcTimeBlocks strataData
                         timeBlockIterationInfo
                         newOverallData
                         newInformationMatrix

computeBreslow
  :: IterationInfo
  -> OverallData
  -> OverallData
  -> Matrix Double
  -> (OverallData, Matrix Double)
computeBreslow iterationInfo overallData tiedData informationMatrix
  | iterationInfo.nEvents == 0 =
    let newLogLikelihood = overallData.logLikelihood + tiedData.logLikelihood
        newOverallData = (overallData :: OverallData) { logLikelihood = newLogLikelihood }
    in  (newOverallData, informationMatrix)
  | otherwise =
    let
        newSumWeightedRisk = overallData.sumWeightedRisk
                             + tiedData.sumWeightedRisk
        newXBarUnscaled = add overallData.xBarUnscaled
                              tiedData.xBarUnscaled
        newXBar = scale (1 / newSumWeightedRisk) newXBarUnscaled
        newOverallData = OverallData
          { sumWeights = 0
          , sumWeightedRisk = newSumWeightedRisk
          , logLikelihood = overallData.logLikelihood
                            + tiedData.logLikelihood
                            - (tiedData.sumWeights * log newSumWeightedRisk)
          , score = add overallData.score (scale (- tiedData.sumWeights) newXBar)
          , xBarUnscaled = newXBarUnscaled
          , informationTerm1 = add overallData.informationTerm1
                                   tiedData.informationTerm1
          }
        blockInformation = scale (tiedData.sumWeights / newSumWeightedRisk)
                                 (add newOverallData.informationTerm1
                                      (scale (- newSumWeightedRisk)
                                             (outer newXBar newXBar)))
        newInformationMatrix = add informationMatrix blockInformation
    in  (newOverallData, newInformationMatrix)

computeEfron
  :: IterationInfo
  -> OverallData
  -> OverallData
  -> Matrix Double
  -> (OverallData, Matrix Double)
computeEfron iterationInfo overallData tiedData informationMatrix
  | iterationInfo.nEvents <= 1 =
    computeBreslow iterationInfo overallData tiedData informationMatrix
  | otherwise =
    undefined -- FIXME

calcTimeBlocksSubjects
  :: StrataData
  -> IterationInfo
  -> OverallData
  -> OverallData
  -> (IterationInfo, OverallData, OverallData)
calcTimeBlocksSubjects strataData iterationInfo overallData tiedData
  -- Case: we've either seen all of the subjects or we've found a subject with a
  -- different censoring or event time. This is the base case
  | (iterationInfo.subjectIndex < 0) -- TODO: create helper functions for checks
      || checkDifferentStrata strataData iterationInfo
      || (strataData.time ! iterationInfo.subjectIndex) /= iterationInfo.time =
      (iterationInfo, overallData, tiedData)
  -- Case: the current subject is part of the current censoring or event tied
  -- time block (note that it could be the first subject in the block). We
  -- calculate and accumulate all of the relevant terms for the subject and
  -- recursively call `calcTimeBlocksSubjects` to conditionally compute the next
  -- subject in the time block
  | otherwise =
      let subjectWeight = strataData.weights ! iterationInfo.subjectIndex
          subjectXProdBeta = strataData.xProdBeta ! iterationInfo.subjectIndex
          subjectWeightedRisk = subjectWeight * exp subjectXProdBeta
          subjectXRow = strataData.xDesignMatrix ! iterationInfo.subjectIndex
          weightedSubjectXRow = scale subjectWeightedRisk subjectXRow
          newInformationTerm1 = scale subjectWeightedRisk
                                      (outer subjectXRow subjectXRow)
      in  case strataData.eventStatus V.! iterationInfo.subjectIndex of
            Censored ->
              let newOverallData = overallData
                    { sumWeightedRisk = overallData.sumWeightedRisk
                                        + subjectWeightedRisk
                    , xBarUnscaled = add overallData.xBarUnscaled
                                         weightedSubjectXRow
                    , informationTerm1 = newInformationTerm1
                    }
                  newIterationInfo = iterationInfo
                    { subjectIndex = iterationInfo.subjectIndex - 1
                    }
              in  calcTimeBlocksSubjects
                    strataData
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
                    , score = add tiedData.score weightedSubjectXRow
                    , xBarUnscaled = add tiedData.xBarUnscaled
                                         weightedSubjectXRow
                    , informationTerm1 = newInformationTerm1
                    }
                  newIterationInfo = iterationInfo
                    { nEvents = iterationInfo.nEvents + 1
                    , subjectIndex = iterationInfo.subjectIndex - 1
                    }
              in  calcTimeBlocksSubjects
                    strataData
                    newIterationInfo
                    overallData
                    newTiedData

-- addSubject :: VS.Vector Double -> Int -> ()


-- need to update:
--   weightedRisks
--   logLikelihood
--   score
--   informationTerm1

-- testMatrix :: Matrix Double
-- testMatrix = matrix 2 [1..4]

v1 :: Vector Double
v1 = L.fromList [3, 6]

v2 :: Vector Double
v2 = L.fromList [3, 6]

m :: Matrix Double
m = fromColumns [v1, v2]

-- testResult :: Vector Double
-- testResult =  testMatrix #> testVector

createInitialData :: Int -> OverallData
createInitialData p
  = OverallData
      { sumWeights = 0
      , sumWeightedRisk = 0
      , logLikelihood = 0
      , score = VS.replicate p 0
      , xBarUnscaled = VS.replicate p 0
      , informationTerm1 = createEmptyMatrix p
      }

createEmptyMatrix :: Int -> Matrix Double
createEmptyMatrix p = diag (VS.replicate p 0)

checkDifferentStrata :: StrataData -> IterationInfo -> Bool
checkDifferentStrata strataData iterationInfo =
  let subjectStratum = strataData.stratum VS.! iterationInfo.subjectIndex
  in  subjectStratum /= iterationInfo.stratum
