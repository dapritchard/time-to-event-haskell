{-# LANGUAGE DuplicateRecordFields #-}
-- |

module CoxPH.CoxPHUpdateNewtonRaphson where

import CoxPH.Data
-- import Numeric.LinearAlgebra.Data -- (Matrix, Vector, fromColumns)
import Numeric.LinearAlgebra qualified as L -- (Matrix, Vector, fromColumns)
import Numeric.LinearAlgebra ( (!), (#>), Matrix, Vector, add, outer, scale, fromColumns )
import Data.Vector.Storable qualified as VS
import Data.Vector.Unboxed qualified as VU
import Data.Vector qualified as V
import qualified CoxPH.CoxPHUpdateNewtonRaphson as overallData

data StrataData = StrataData {
    time :: VS.Vector Double
  , eventStatus :: V.Vector Delta
  , xDesignMatrix :: Matrix Double
  , xOffset :: Vector Double
  }

updateNewtonRaphson :: [StrataData] -> VS.Vector Double -> Vector Double -> ()
updateNewtonRaphson strataDatas weights beta =
  -- let z = add ()
  ()

updateStrata :: StrataData -> VS.Vector Double -> Vector Double -> ()
updateStrata strataData weights beta =
  let xProdBeta = add (strataData.xDesignMatrix #> beta) strataData.xOffset
      -- risk =
  in  ()
-- updateTime

-- data TTEData = TTEData {

--   }

data IterationInfo = IterationInfo {
    nSubjects :: Int
  , subjectIndex :: Int
  , nDead :: Int
  , time :: Double
  }

data OverallData = OverallData {
    sumWeightedRisk :: Double
  , logLikelihood :: Double
  , score :: Vector Double
  , xBarUnscaled :: Vector Double
  , informationUnscaled :: Matrix Double
  }

data TiedData = TiedData {
    sumWeights :: Double
  , sumWeightedRisk :: Double
  , logLikelihood :: Double
  , score :: Vector Double
  , xBarUnscaled :: Vector Double
  , informationUnscaled :: Matrix Double
  }

calcTimeBlock :: Int -> Int -> VS.Vector Double -> ()
calcTimeBlock nSubjects subjectIndex xProdBeta
  | subjectIndex >= nSubjects = ()
  | otherwise =
    ()

-- calcSubject :: Int -> Int -> Int -> Double -> Vector Double -> StrataData -> OverallData -> TiedData -> ()
-- calcSubject nSubjects subjectIndex nDead time weightedRisk strataData overallData tiedData
calcSubject :: IterationInfo -> Vector Double -> StrataData -> OverallData -> TiedData
               -> (IterationInfo, StrataData, OverallData, TiedData)
calcSubject iterationInfo weightedRisk strataData overallData tiedData
  -- case: we've either seen all of the subjects or we've found a subject with a
  -- different censoring or event time. This is the base case
  | (iterationInfo.subjectIndex >= iterationInfo.nSubjects)
    || ((strataData.time ! iterationInfo.subjectIndex) /= iterationInfo.time) =
      (iterationInfo, strataData, overallData, tiedData)
  -- case: TODO
  | otherwise =
      let subjectWeightedRisk = weightedRisk ! iterationInfo.subjectIndex
          subjectXRow = strataData.xDesignMatrix ! iterationInfo.subjectIndex
          weightedSubjectXRow = scale subjectWeightedRisk subjectXRow
          subjectCovariateTerm = scale subjectWeightedRisk
                                       (outer subjectXRow
                                              subjectXRow)
      in  case strataData.eventStatus V.! iterationInfo.subjectIndex of
            Censored ->
              let newSumWeightedRisk = (overallData.sumWeightedRisk
                                        + subjectWeightedRisk)
                  newXBarUnscaled = add overallData.xBarUnscaled
                                        weightedSubjectXRow
                  newOverallData = OverallData {
                        sumWeightedRisk = newSumWeightedRisk
                      , logLikelihood = overallData.logLikelihood
                      , score = overallData.score
                      , xBarUnscaled = newXBarUnscaled
                      , informationUnscaled = overallData.informationUnscaled
                    }
              in  (iterationInfo, strataData, newOverallData, tiedData)
            ObservedEvent ->
              let newSumWeights = tiedData.sumWeights + subjectWeightedRisk
                  newSumWeightedRisk = tiedData.sumWeights + subjectWeightedRisk
                  newLogLikelihood = tiedData.logLikelihood + subjectWeightedRisk
                  newScore = add tiedData.score subjectXRow
                  newXBarUnscaled = add tiedData.xBarUnscaled weightedSubjectXRow
                  newTiedData = tiedData {
                          sumWeights = newSumWeights
                        , sumWeightedRisk = newSumWeightedRisk
                        , logLikelihood = newLogLikelihood
                        , score = newScore
                        , xBarUnscaled = newXBarUnscaled
                        , informationUnscaled = tiedData.informationUnscaled
                        }
              in  (iterationInfo, strataData, overallData, newTiedData)

-- addSubject :: VS.Vector Double -> Int -> ()


-- need to update:
--   weightedRisk
--   logLikelihood
--   score
--   informationUnscaled

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
