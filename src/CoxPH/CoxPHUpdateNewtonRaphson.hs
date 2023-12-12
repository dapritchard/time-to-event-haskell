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
  , nDead :: Int
  , time :: Double
  }

data StrataData = StrataData {
    time :: VS.Vector Double
  , eventStatus :: V.Vector Delta
  , xDesignMatrix :: Matrix Double
  , xOffset :: Vector Double
  , weights :: Vector Double
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
      weightedRisk = VS.zipWith (*) strataData.weights xProdBeta
  in  calcTimeBlock strataData iterationInfo weightedRisk overallData

calcTimeBlock
  :: StrataData
  -> IterationInfo
  -> Vector Double
  -> OverallData
  -> (IterationInfo, OverallData)
calcTimeBlock strataData iterationInfo weightedRisk overallData =
  let p = VS.length overallData.score
      initialTiedData = createInitialTiedData p
      (newIterationInfo, newOverallData, newTiedData) =
        calcTimeBlockSubjects
          strataData iterationInfo weightedRisk overallData initialTiedData
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
  -> Vector Double
  -> OverallData
  -> TiedData
  -> (IterationInfo, OverallData, TiedData)
calcTimeBlockSubjects strataData iterationInfo weightedRisk overallData tiedData
  -- Case: we've either seen all of the subjects or we've found a subject with a
  -- different censoring or event time. This is the base case
  | (iterationInfo.subjectIndex >= iterationInfo.nSubjects)
      || (strataData.time ! iterationInfo.subjectIndex) /= iterationInfo.time =
      (iterationInfo, overallData, tiedData)
  -- Case: the current subject is part of the censoring or event tied time block
  -- (note that it could be the first subject in the block)
  | otherwise =
      let subjectWeightedRisk = weightedRisk ! iterationInfo.subjectIndex
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
                  newIterationInfo = iterationInfo {
                      subjectIndex = iterationInfo.subjectIndex + 1
                    }
              in  calcTimeBlockSubjects
                    strataData
                    newIterationInfo
                    weightedRisk
                    newOverallData
                    tiedData
            ObservedEvent ->
              let newTiedData = TiedData {
                      sumWeights = tiedData.sumWeights + subjectWeightedRisk -- FIXME: wrong!!
                    , sumWeightedRisk = tiedData.sumWeights
                                        + subjectWeightedRisk
                    , logLikelihood = tiedData.logLikelihood
                                      + subjectWeightedRisk
                    , score = add tiedData.score subjectXRow
                    , xBarUnscaled = add tiedData.xBarUnscaled
                                         weightedSubjectXRow
                    , informationTerm1 = newInformationTerm1
                    }
                  newIterationInfo = iterationInfo {
                      nDead = iterationInfo.nDead + 1
                    , subjectIndex = iterationInfo.subjectIndex + 1
                    }
              in  calcTimeBlockSubjects
                    strataData
                    newIterationInfo
                    weightedRisk
                    overallData
                    newTiedData

-- addSubject :: VS.Vector Double -> Int -> ()


-- need to update:
--   weightedRisk
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
