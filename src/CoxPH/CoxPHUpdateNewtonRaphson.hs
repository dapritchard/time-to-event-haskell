{-# LANGUAGE DuplicateRecordFields #-}
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

data IterationInfo = IterationInfo {
    nSubjects :: Int
  , subjectIndex :: Int
  , nEvents :: Int
  , time :: Double
  }

data StrataData = StrataData {
    time :: VS.Vector Double
  , eventStatus :: V.Vector Delta
  , xDesignMatrix :: Matrix Double
  , xOffset :: Vector Double
  , weights :: Vector Double
  , xProdBeta :: Vector Double
  }

data OverallData = OverallData {
    sumWeightedRisk :: Double
  , logLikelihood :: Double
  , score :: Vector Double
  , xBarUnscaled :: Vector Double
  , informationTerm1 :: Matrix Double
  }

data TiedData = TiedData {
    sumWeights :: Double
  , sumWeightedRisk :: Double
  , logLikelihood :: Double
  , score :: Vector Double
  , xBarUnscaled :: Vector Double
  , informationTerm1 :: Matrix Double
  }

updateStrata
  :: Vector Double
  -> StrataData
  -> IterationInfo
  -> OverallData
  -> (IterationInfo, OverallData)
updateStrata beta strataData iterationInfo overallData =
  let xProdBeta = add (strataData.xDesignMatrix #> beta) strataData.xOffset
      weightedRisks = VS.zipWith calcWeightedRisk strataData.weights xProdBeta
  in  calcTimeBlock strataData iterationInfo overallData
  where
    calcWeightedRisk :: Double -> Double -> Double
    calcWeightedRisk w x = w * exp x

calcTimeBlock
  :: StrataData
  -> IterationInfo
  -> OverallData
  -> (IterationInfo, OverallData)
calcTimeBlock strataData iterationInfo overallData =
  let p = VS.length overallData.score
      initialTiedData = createInitialTiedData p
      (newIterationInfo, newOverallData, newTiedData) =
        calcTimeBlockSubjects
          strataData iterationInfo overallData initialTiedData
      newXBarUnscaled = newOverallData.xBarUnscaled + newTiedData.xBarUnscaled
      newXBar = scale newTiedData.sumWeights newXBarUnscaled
      updatedOverallData = OverallData {
          sumWeightedRisk = newOverallData.sumWeightedRisk
                            + newTiedData.sumWeightedRisk
        , logLikelihood = newOverallData.logLikelihood
                          - (newTiedData.sumWeights
                             * log newTiedData.logLikelihood)
        , score = add newOverallData.score
                      (scale newTiedData.sumWeights newXBar)
        , xBarUnscaled = newXBarUnscaled
        , informationTerm1 = add newOverallData.informationTerm1
                                 newTiedData.informationTerm1
        }
  in  (newIterationInfo, updatedOverallData) -- FIXME: need to update information

calcTimeBlockSubjects
  :: StrataData
  -> IterationInfo
  -> OverallData
  -> TiedData
  -> (IterationInfo, OverallData, TiedData)
calcTimeBlockSubjects strataData iterationInfo overallData tiedData
  -- Case: we've either seen all of the subjects or we've found a subject with a
  -- different censoring or event time. This is the base case
  | (iterationInfo.subjectIndex < 0)                             -- FIXME: need to check that we don't cross strata
      || (strataData.time ! iterationInfo.subjectIndex) /= iterationInfo.time =
      (iterationInfo, overallData, tiedData)
  -- Case: the current subject is part of the censoring or event tied time block
  -- (note that it could be the first subject in the block). We calculate and
  -- accumulate all of the relevant terms for the subject and recursively call
  -- `calcTimeBlockSubjects` to conditionally compute the next subject in the
  -- time block
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
              let newOverallData = OverallData
                    { sumWeightedRisk = overallData.sumWeightedRisk
                                        + subjectWeightedRisk
                    , logLikelihood = overallData.logLikelihood
                    , score = overallData.score
                    , xBarUnscaled = add overallData.xBarUnscaled
                                         weightedSubjectXRow
                    , informationTerm1 = newInformationTerm1
                    }
                  newIterationInfo = iterationInfo
                    { subjectIndex = iterationInfo.subjectIndex - 1
                    }
              in  calcTimeBlockSubjects
                    strataData
                    newIterationInfo
                    newOverallData
                    tiedData
            ObservedEvent ->
              let newTiedData = TiedData
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
              in  calcTimeBlockSubjects
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

createInitialTiedData :: Int -> TiedData
createInitialTiedData p
  = TiedData {
        sumWeights = 0
      , sumWeightedRisk = 0
      , logLikelihood = 0
      , score = VS.replicate p 0
      , xBarUnscaled = VS.replicate p 0
      , informationTerm1 = createEmptyMatrix p
      }

createEmptyMatrix :: Int -> Matrix Double
createEmptyMatrix p = diag (VS.replicate p 0)
