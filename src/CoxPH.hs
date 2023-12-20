{-# LANGUAGE MultiWayIf #-}

module CoxPH (
  CoxPHConvergenceFailure(..),
  CoxPHMethod(..),
  Delta(..),
  LastInStrataIndicator(..),
  ScaleCovariateIndicator(..),
  coxph
  ) where

import CoxPH.CenterAndScale (centerAndScaleCovs)
import CoxPH.CoxPHUpdateNewtonRaphson ( NRResults(..), TTEData(..), coxPHUpdateNewtonRaphson )
import CoxPH.Data
import Data.Vector.Storable qualified as VS
import Data.Vector qualified as V
import Data.Text qualified as T
import Numeric.LinearAlgebra qualified as L

-- type CoxPHResult = Either CoxPHConvergenceFailure (VU.Vector Double)
type CoxPHResult = Either T.Text (VS.Vector Double)

coxph :: VS.Vector Double                  -- length n, the per-subject minumum times of the event or censoring times
      -> V.Vector Delta                    -- length n, indicators for whether patient observered an event or censoring
      -> V.Vector (VS.Vector Double)       -- dimension n by p, X design matrix
      -> VS.Vector Double                  -- length n, X offset vector
      -> VS.Vector Double                  -- length n, subject weights
      -> VS.Vector Int                     -- length n, strata
      -> VS.Vector Double                  -- length p, starting values for the beta coefficients
      -> V.Vector ScaleCovariateIndicator  -- length p, indicators for whether a given covariate should be centered and scaled
      -> CoxPHMethod                       -- which of the Breslow or Effron methods to use in the event of tied event times
      -> Int                               -- the maximumn number of Newton-Rhapson iterations to perform
      -> Double                            -- the value for which the absolute value of 1 minus ratio of the likelihood for two consecutive iterations must be below for convergence to be achieved
      -- -> Double                            -- TODO
      -> CoxPHResult
coxph times
      eventStatuses
      xDesignDataFrame
      xOffset
      weights
      strata
      beta
      scaleIndicators
      tiesMethod
      maxIterations
      epsilon
      -- tolerance
      = do
  -- in Right (VS.fromList [])
  centeredAndScaledCovsResults <- centerAndScaleCovs xDesignDataFrame
                                                     weights
                                                     scaleIndicators
  let (centeredAndScaledCovs, scales) = V.unzip centeredAndScaledCovsResults
      -- scales = V.map snd centeredAndScaledCovsResults
      xDesignMatrix = L.fromColumns (V.toList centeredAndScaledCovs)
      tteData = TTEData
        { time = times
        , eventStatus = eventStatuses
        , xDesignMatrix = xDesignMatrix
        , stratum = strata
        , weights = weights
        , xProdBeta = L.add xOffset (xDesignMatrix L.#> beta)
        , tiesMethod = tiesMethod
        }
  (betaOffset, nrResults) <- calcBetaOffset tteData
  let newBeta = L.add beta betaOffset
  finalBeta <- updateStep maxIterations
                     epsilon
                     2 -- since we've manually done one iteration so far
                     newBeta
                     xOffset
                     nrResults.sumLogLikelihood
                     (updateTTEData newBeta xOffset tteData)
  Right (VS.convert scales / finalBeta)

updateStep :: Int -> Double -> Int -> VS.Vector Double -> VS.Vector Double -> Double -> TTEData -> Either T.Text (VS.Vector Double)
updateStep maxIterations epsilon iteration beta xOffset logLikelihood tteData
  | iteration > maxIterations =
    Left "Did not converge in the alloted number of iterations"
  | otherwise =
    case calcBetaOffset tteData of
      Left e -> Left e
      Right (betaOffset, nrResults) ->
        let newBeta = L.add beta betaOffset
        in  if -- Case: we've achieved convergence;
               | checkNRConvergence epsilon logLikelihood nrResults ->
                 Right newBeta
               -- Case: the logLikelihood is moving in the wrong direction
               | checkDecreasingLogLikelihood logLikelihood nrResults ->
                 Left "Likelihood not nondecreasing"
               -- Case: we haven't achieved convergence yet, but the log
               -- likelihood is moving in the right direction
               | otherwise ->
                 updateStep maxIterations
                            epsilon
                            (iteration + 1)
                            newBeta
                            xOffset
                            nrResults.sumLogLikelihood
                            (updateTTEData newBeta xOffset tteData)

updateTTEData :: VS.Vector Double -> VS.Vector Double -> TTEData -> TTEData
updateTTEData beta xOffset tteData =
  let newXProdBeta = L.add xOffset
                           (tteData.xDesignMatrix L.#> beta)
  in  tteData { xProdBeta = newXProdBeta }

-- calculates I^{-1}(\hat{\beta}^{(n)})\, U(\hat{\beta}^{(n)}), i.e. the term
-- that we add to \hat{\beta}^{(n)} to update it
calcBetaOffset :: TTEData -> Either T.Text (VS.Vector Double, NRResults)
calcBetaOffset tteData = do
  let nrResults = coxPHUpdateNewtonRaphson tteData
      inverseInformationMatrix = L.inv nrResults.informationMatrix  -- FIXME: this can fail. Need to catch failure as appropriate
      betaOffset = inverseInformationMatrix L.#> nrResults.score
  Right (betaOffset, nrResults)

checkNRConvergence :: Double -> Double -> NRResults -> Bool
checkNRConvergence epsilon logLikelihood nrResults =
  abs (1 - (logLikelihood / nrResults.sumLogLikelihood)) < epsilon

checkDecreasingLogLikelihood :: Double -> NRResults -> Bool
checkDecreasingLogLikelihood logLikelihood nrResults =
  logLikelihood >= nrResults.sumLogLikelihood

result :: CoxPHResult
result = coxph (VS.fromList [1, 1, 2, 2, 3, 3, 4]) -- times
               (V.fromList [ObservedEvent, Censored, ObservedEvent, ObservedEvent, ObservedEvent, Censored, ObservedEvent])
               (V.fromList [ VS.fromList [1, 1, 1, 0, 2, 0, 0]
                           , VS.fromList [0, 0, 1, 1, 0, 1, 0]
                           ]) -- xDesignDataFrame
               (VS.fromList [0, 0, 0, 0, 0, 0, 0]) -- xOffset
               (VS.fromList [1, 1, 1, 1, 1, 1, 1]) -- weights
               (VS.fromList [0, 0, 0, 0, 0, 0, 0]) -- strata
               (VS.fromList [0, 0]) -- beta
               (V.fromList [ScaleCovariateYes, ScaleCovariateNo]) -- scaleIndicators
               Breslow -- tiesMethod
               20 -- maxIterations
               0.000000000001818989 -- epsilon

-- data Delta = ObservedEvent | Censored

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
