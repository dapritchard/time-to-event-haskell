{-# LANGUAGE DuplicateRecordFields #-}

module TestCoxPH.TestCoxPHUpdateNewtonRaphson where

import CoxPH.CenterAndScale
import CoxPH.CoxPHUpdateNewtonRaphson
import CoxPH.Data
import Data.Vector qualified as V
import Data.Vector.Storable qualified as VS
import Data.Either (fromRight)
import Numeric.LinearAlgebra qualified as L

-- Create test input data ------------------------------------------------------

-- The following are the data inputs resulting from the following call
-- immediately before being input to `Ccoxfit6_iter`
--
-- test1 <- tibble(
--   time   = c(1, 2, 3, 4, 5, 6, 7),
--   status = c(1, 0, 1, 1, 1, 0, 1),
--   x      = c(1, 1, 1, 0, 2, 0, 0),
--   sex    = c(0, 0, 1, 1, 0, 1, 0)
-- )
-- -- test2 <- list(time=c(4,3,1,1,2,2,3),
-- --                status=c(1,1,1,0,1,1,0),
-- --                x=c(0,2,1,1,1,0,0),
-- --                sex=c(0,0,0,0,1,1,1))
-- coxph(Surv(time, status) ~ x + sex, test1)
--
-- Browse[2]> > as.integer(maxiter)
-- [1] 20
-- Browse[2]> stime
-- [1] 1 1 2 2 3 3 4
-- Browse[2]> sstat
-- [1] 1 0 1 1 1 0 1
-- Browse[2]> > x[sorted,]
--   x sex
-- 3 1   0
-- 4 1   0
-- 5 1   1
-- 6 0   1
-- 2 2   0
-- 7 0   1
-- 1 0   0
-- Browse[2]> > as.double(offset[sorted])
-- [1] 0 0 0 0 0 0 0
-- Browse[2]> weights
-- [1] 1 1 1 1 1 1 1
-- Browse[2]> newstrat
-- [1] 0 0 0 0 0 0 0
-- Browse[2]> > as.integer(method=="efron")
-- [1] 0
-- Browse[2]> > as.double(control$eps)
-- [1] 0.000000001
-- Browse[2]> > as.double(control$toler.chol)
-- [1] 0.000000000001818989
-- Browse[2]> > as.vector(init)
-- [1] 0 0
-- Browse[2]> > ifelse(zero.one, 0L, 1L)
--   x sex
--   1   0

times :: VS.Vector Double
times = VS.fromList [0, 1, 2, 3, 4, 5, 6]
eventStatuses :: V.Vector Delta
eventStatuses = V.fromList [ ObservedEvent
                           , Censored
                           , ObservedEvent
                           , ObservedEvent
                           , ObservedEvent
                           , Censored
                           , ObservedEvent
                           ]
xDesignDataFrame :: V.Vector (VS.Vector Double)
xDesignDataFrame = V.fromList [ VS.fromList [1, 1, 1, 0, 2, 0, 0]
                              , VS.fromList [0, 0, 1, 1, 0, 1, 0]
                              ]
xOffset :: VS.Vector Double
xOffset = VS.fromList [0, 0, 0, 0, 0, 0, 0]
testWeights :: VS.Vector Double
testWeights = VS.fromList [1, 1, 1, 1, 1, 1, 1]
strata :: VS.Vector Int
strata = VS.fromList [0, 0, 0, 0, 0, 0, 0]
initialBeta :: VS.Vector Double
initialBeta = VS.fromList [0, 0]
scaleIndicators :: V.Vector ScaleCovariateIndicator
-- scaleIndicators = V.fromList [ScaleCovariateYes, ScaleCovariateNo]
scaleIndicators = V.fromList [ScaleCovariateNo, ScaleCovariateNo]
testTiesMethod :: CoxPHMethod
testTiesMethod = Breslow
maxIterations :: Int
maxIterations = 20
epsilon :: Double
epsilon = 0.000000000001818989

-- times :: VS.Vector Double
-- times = VS.fromList [1, 1, 2, 2, 3, 3, 4]
-- eventStatuses :: V.Vector Delta
-- eventStatuses = V.fromList [ ObservedEvent
--                            , Censored
--                            , ObservedEvent
--                            , ObservedEvent
--                            , ObservedEvent
--                            , Censored
--                            , ObservedEvent
--                            ]
-- xDesignDataFrame :: V.Vector (VS.Vector Double)
-- xDesignDataFrame = V.fromList [ VS.fromList [1, 1, 1, 0, 2, 0, 0]
--                               , VS.fromList [0, 0, 1, 1, 0, 1, 0]
--                               ]
-- xOffset :: VS.Vector Double
-- xOffset = VS.fromList [0, 0, 0, 0, 0, 0, 0]
-- testWeights :: VS.Vector Double
-- testWeights = VS.fromList [1, 1, 1, 1, 1, 1, 1]
-- strata :: VS.Vector Int
-- strata = VS.fromList [0, 0, 0, 0, 0, 0, 0]
-- initialBeta :: VS.Vector Double
-- initialBeta = VS.fromList [0, 0]
-- scaleIndicators :: V.Vector ScaleCovariateIndicator
-- scaleIndicators = V.fromList [ScaleCovariateYes, ScaleCovariateNo]
-- testTiesMethod :: CoxPHMethod
-- testTiesMethod = Breslow
-- maxIterations :: Int
-- maxIterations = 20
-- epsilon :: Double
-- epsilon = 0.000000000001818989

centeredAndScaledCovsResults :: V.Vector (VS.Vector Double, Double)
centeredAndScaledCovsResults =
  fromRight (V.singleton (VS.singleton 0, 0))
            (centerAndScaleCovs xDesignDataFrame
                                testWeights
                                scaleIndicators)

centeredAndScaledCovs :: V.Vector (VS.Vector Double)
(centeredAndScaledCovs, _) = V.unzip centeredAndScaledCovsResults

testXDesignMatrix :: L.Matrix Double
testXDesignMatrix = L.fromColumns (V.toList xDesignDataFrame)

tteData :: TTEData
tteData = TTEData
  { time = times
  , eventStatus = eventStatuses
  , xDesignMatrix = testXDesignMatrix
  , stratum = strata
  , weights = testWeights
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta)
  , tiesMethod = testTiesMethod
  }

nrResults :: NRResults
nrResults = coxPHUpdateNewtonRaphson tteData

-- -- Patient 6
-- (IterationInfo {subjectIndex = 5, time = 7.0, stratum = 0, nEvents = 1},NRTerms {sumWeights = 0.0, sumWeightedRisk = 0.0, logLikelihood = 0.0, score = [0.0,0.0], xBarUnscaled = [0.0,0.0], informationTerm1 = (2><2)
--  [ 0.0, 0.0
--  , 0.0, 0.0 ]},NRTerms {sumWeights = 1.0, sumWeightedRisk = 1.0, logLikelihood = 0.0, score = [0.0,0.0], xBarUnscaled = [0.0,0.0], informationTerm1 = (2><2)
--  [ 0.0, 0.0
--  , 0.0, 0.0 ]})

tteData6 :: TTEData
tteData6 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [6]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [6]
  , xDesignMatrix = testXDesignMatrix L.? [6]
  , stratum = strata `VS.backpermute` VS.fromList [6]
  , weights = testWeights `VS.backpermute` VS.fromList [6]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [6]
  , tiesMethod = testTiesMethod
  }

out6 :: (IterationInfo, NRTerms, L.Matrix Double)
out6 = calcTimeBlocks tteData6
                      (IterationInfo 0 0 0 0)
                      (createInitialData 2)
                      (createEmptyMatrix 2)

tteData5 :: TTEData
tteData5 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [5]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [5]
  , xDesignMatrix = testXDesignMatrix L.? [5]
  , stratum = strata `VS.backpermute` VS.fromList [5]
  , weights = testWeights `VS.backpermute` VS.fromList [5]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [5]
  , tiesMethod = testTiesMethod
  }

out5 :: (IterationInfo, NRTerms, L.Matrix Double)
out5 = calcTimeBlocks tteData5
                      ((fst3 out6) { subjectIndex = 0 })
                      (snd3 out6)
                      (thd3 out6)

tteData4 :: TTEData
tteData4 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [4]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [4]
  , xDesignMatrix = testXDesignMatrix L.? [4]
  , stratum = strata `VS.backpermute` VS.fromList [4]
  , weights = testWeights `VS.backpermute` VS.fromList [4]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [4]
  , tiesMethod = testTiesMethod
  }

out4 :: (IterationInfo, NRTerms, L.Matrix Double)
out4 = calcTimeBlocks tteData4
                      ((fst3 out5) { subjectIndex = 0 })
                      (snd3 out5)
                      (thd3 out5)

uncurry3 :: (a -> b -> c -> d) -> (a, b, c) -> d
uncurry3 f (a, b, c) = f a b c

-- z :: (IterationInfo, NRTerms, NRTerms)
-- z = calcTimeBlocksSubjects tteData
--                            (IterationInfo 5 5 0 1)
--                            (snd3 out6)
--                            (createInitialData 2)

fst3 :: (a, b, c) -> a
fst3 (a, _, _) = a

snd3 :: (a, b, c) -> b
snd3 (_, b, _) = b

thd3 :: (a, b, c) -> c
thd3 (_, _, c) = c
