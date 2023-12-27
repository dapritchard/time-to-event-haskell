{-# LANGUAGE DuplicateRecordFields #-}

module TestCoxPH.TestCoxPHUpdateNewtonRaphson where

import CoxPH
import CoxPH.CenterAndScale
import CoxPH.CoxPHUpdateNewtonRaphson
import CoxPH.Data
import Data.Text qualified as T
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
scaleIndicators = V.fromList [ScaleCovariateYes, ScaleCovariateNo]
-- scaleIndicators = V.fromList [ScaleCovariateNo, ScaleCovariateNo]
testTiesMethod :: CoxPHMethod
testTiesMethod = Breslow
maxIterations :: Int
maxIterations = 20
epsilon :: Double
epsilon = 0.000000000001818989

out :: Either T.Text (VS.Vector Double)
out = coxph times
            eventStatuses
            xDesignDataFrame
            xOffset
            testWeights
            strata
            initialBeta
            scaleIndicators
            testTiesMethod
            maxIterations
            epsilon

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

tteData0 :: TTEData
tteData0 = tteData
nrResults0 :: NRResults
nrResults0 = coxPHUpdateNewtonRaphson tteData0
inverseInformationMatrix0 :: L.Matrix Double
inverseInformationMatrix0 = L.inv nrResults0.informationMatrix
beta1 :: L.Vector Double
beta1 = inverseInformationMatrix0 L.#> nrResults0.score

tteData1 :: TTEData
tteData1 = tteData0 { xProdBeta = tteData0.xDesignMatrix L.#> beta1 }
nrResults1 :: NRResults
nrResults1 = coxPHUpdateNewtonRaphson tteData1
inverseInformationMatrix1 :: L.Matrix Double
inverseInformationMatrix1 = L.inv nrResults1.informationMatrix
beta2 :: L.Vector Double
beta2 = L.add beta1
              (inverseInformationMatrix1 L.#> nrResults1.score)

tteData2 :: TTEData
tteData2 = tteData1 { xProdBeta = tteData1.xDesignMatrix L.#> beta2 }
nrResults2 :: NRResults
nrResults2 = coxPHUpdateNewtonRaphson tteData2
inverseInformationMatrix2 :: L.Matrix Double
inverseInformationMatrix2 = L.inv nrResults2.informationMatrix
beta3 :: L.Vector Double
beta3 = L.add beta2
              (inverseInformationMatrix2 L.#> nrResults2.score)

tteData3 :: TTEData
tteData3 = tteData2 { xProdBeta = tteData2.xDesignMatrix L.#> beta3 }
nrResults3 :: NRResults
nrResults3 = coxPHUpdateNewtonRaphson tteData3
inverseInformationMatrix3 :: L.Matrix Double
inverseInformationMatrix3 = L.inv nrResults3.informationMatrix
beta4 :: L.Vector Double
beta4 = L.add beta3
              (inverseInformationMatrix3 L.#> nrResults3.score)

-- -- Patient 6
-- (IterationInfo {subjectIndex = 5, time = 7.0, stratum = 0, nEvents = 1},NRTerms {sumWeights = 0.0, sumWeightedRisk = 0.0, logLikelihood = 0.0, score = [0.0,0.0], xBarUnscaled = [0.0,0.0], informationTerm1 = (2><2)
--  [ 0.0, 0.0
--  , 0.0, 0.0 ]},NRTerms {sumWeights = 1.0, sumWeightedRisk = 1.0, logLikelihood = 0.0, score = [0.0,0.0], xBarUnscaled = [0.0,0.0], informationTerm1 = (2><2)
--  [ 0.0, 0.0
--  , 0.0, 0.0 ]})

tteDataT6 :: TTEData
tteDataT6 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [6]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [6]
  , xDesignMatrix = testXDesignMatrix L.? [6]
  , stratum = strata `VS.backpermute` VS.fromList [6]
  , weights = testWeights `VS.backpermute` VS.fromList [6]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [6]
  , tiesMethod = testTiesMethod
  }

out6 :: (IterationInfo, NRTerms, L.Matrix Double)
out6 = calcTimeBlocks tteDataT6
                      (IterationInfo 0 0 0 0)
                      (createInitialData 2)
                      (createEmptyMatrix 2)

tteDataT5 :: TTEData
tteDataT5 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [5]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [5]
  , xDesignMatrix = testXDesignMatrix L.? [5]
  , stratum = strata `VS.backpermute` VS.fromList [5]
  , weights = testWeights `VS.backpermute` VS.fromList [5]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [5]
  , tiesMethod = testTiesMethod
  }

out5 :: (IterationInfo, NRTerms, L.Matrix Double)
out5 = calcTimeBlocks tteDataT5
                      ((fst3 out6) { subjectIndex = 0 })
                      (snd3 out6)
                      (thd3 out6)

tteDataT4 :: TTEData
tteDataT4 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [4]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [4]
  , xDesignMatrix = testXDesignMatrix L.? [4]
  , stratum = strata `VS.backpermute` VS.fromList [4]
  , weights = testWeights `VS.backpermute` VS.fromList [4]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [4]
  , tiesMethod = testTiesMethod
  }

out4 :: (IterationInfo, NRTerms, L.Matrix Double)
out4 = calcTimeBlocks tteDataT4
                      ((fst3 out5) { subjectIndex = 0 })
                      (snd3 out5)
                      (thd3 out5)

tteDataT3 :: TTEData
tteDataT3 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [3]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [3]
  , xDesignMatrix = testXDesignMatrix L.? [3]
  , stratum = strata `VS.backpermute` VS.fromList [3]
  , weights = testWeights `VS.backpermute` VS.fromList [3]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [3]
  , tiesMethod = testTiesMethod
  }

out3 :: (IterationInfo, NRTerms, L.Matrix Double)
out3 = calcTimeBlocks tteDataT3
                      ((fst3 out4) { subjectIndex = 0 })
                      (snd3 out4)
                      (thd3 out4)

tteDataT2 :: TTEData
tteDataT2 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [2]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [2]
  , xDesignMatrix = testXDesignMatrix L.? [2]
  , stratum = strata `VS.backpermute` VS.fromList [2]
  , weights = testWeights `VS.backpermute` VS.fromList [2]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [2]
  , tiesMethod = testTiesMethod
  }

out2 :: (IterationInfo, NRTerms, L.Matrix Double)
out2 = calcTimeBlocks tteDataT2
                      ((fst3 out3) { subjectIndex = 0 })
                      (snd3 out3)
                      (thd3 out3)

tteDataT1 :: TTEData
tteDataT1 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [1]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [1]
  , xDesignMatrix = testXDesignMatrix L.? [1]
  , stratum = strata `VS.backpermute` VS.fromList [1]
  , weights = testWeights `VS.backpermute` VS.fromList [1]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [1]
  , tiesMethod = testTiesMethod
  }

out1 :: (IterationInfo, NRTerms, L.Matrix Double)
out1 = calcTimeBlocks tteDataT1
                      ((fst3 out2) { subjectIndex = 0 })
                      (snd3 out2)
                      (thd3 out2)

tteDataT0 :: TTEData
tteDataT0 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [0]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [0]
  , xDesignMatrix = testXDesignMatrix L.? [0]
  , stratum = strata `VS.backpermute` VS.fromList [0]
  , weights = testWeights `VS.backpermute` VS.fromList [0]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> initialBeta) `VS.backpermute` VS.fromList [0]
  , tiesMethod = testTiesMethod
  }

out0 :: (IterationInfo, NRTerms, L.Matrix Double)
out0 = calcTimeBlocks tteDataT0
                      ((fst3 out1) { subjectIndex = 0 })
                      (snd3 out1)
                      (thd3 out1)

beta11 :: L.Vector Double
beta11 = VS.fromList [1, 1]

tteDataT6B11 :: TTEData
tteDataT6B11 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [6]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [6]
  , xDesignMatrix = testXDesignMatrix L.? [6]
  , stratum = strata `VS.backpermute` VS.fromList [6]
  , weights = testWeights `VS.backpermute` VS.fromList [6]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> beta11) `VS.backpermute` VS.fromList [6]
  , tiesMethod = testTiesMethod
  }

out6B11 :: (IterationInfo, NRTerms, L.Matrix Double)
out6B11 = calcTimeBlocks tteDataT6B11
                         (IterationInfo 0 0 0 0)
                         (createInitialData 2)
                         (createEmptyMatrix 2)

tteDataT5B11 :: TTEData
tteDataT5B11 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [5]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [5]
  , xDesignMatrix = testXDesignMatrix L.? [5]
  , stratum = strata `VS.backpermute` VS.fromList [5]
  , weights = testWeights `VS.backpermute` VS.fromList [5]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> beta11) `VS.backpermute` VS.fromList [5]
  , tiesMethod = testTiesMethod
  }

out5B11 :: (IterationInfo, NRTerms, L.Matrix Double)
out5B11 = calcTimeBlocks tteDataT5B11
                         ((fst3 out6B11) { subjectIndex = 0 })
                         (snd3 out6B11)
                         (thd3 out6B11)

tteDataT4B11 :: TTEData
tteDataT4B11 = TTEData
  { time = tteData.time `VS.backpermute` VS.fromList [4]
  , eventStatus = eventStatuses `V.backpermute` V.fromList [4]
  , xDesignMatrix = testXDesignMatrix L.? [4]
  , stratum = strata `VS.backpermute` VS.fromList [4]
  , weights = testWeights `VS.backpermute` VS.fromList [4]
  , xProdBeta = L.add xOffset (testXDesignMatrix L.#> beta11) `VS.backpermute` VS.fromList [4]
  , tiesMethod = testTiesMethod
  }

out4B11 :: (IterationInfo, NRTerms, L.Matrix Double)
out4B11 = calcTimeBlocks tteDataT4B11
                         ((fst3 out5B11) { subjectIndex = 0 })
                         (snd3 out5B11)
                         (thd3 out5B11)

-- outPass1 :: (IterationInfo, NRTerms, L.Matrix Double)
-- outPass = calcTimeBlocks tteDataT6
--                       (IterationInfo 0 0 0 0)
--                       (createInitialData 2)
--                       (createEmptyMatrix 2)

fst3 :: (a, b, c) -> a
fst3 (a, _, _) = a

snd3 :: (a, b, c) -> b
snd3 (_, b, _) = b

thd3 :: (a, b, c) -> c
thd3 (_, _, c) = c
