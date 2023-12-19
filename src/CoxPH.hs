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
      -> VS.Vector Int                     -- length n, indicators for whether a patient had the latest minumum event of censoring time within a strata
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
      stratums
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
  let centeredAndScaledCovs = V.map fst centeredAndScaledCovsResults
      -- scales = V.map snd centeredAndScaledCovsResults
      xDesignMatrix = L.fromColumns (V.toList centeredAndScaledCovs)
      tteData = TTEData
        { time = times
        , eventStatus = eventStatuses
        , xDesignMatrix = xDesignMatrix
        , stratum = stratums
        , weights = weights
        , xProdBeta = L.add xOffset (xDesignMatrix L.#> beta)
        , tiesMethod = tiesMethod
        }
  (betaOffset, nrResults) <- calcBetaOffset tteData
  let newBeta = L.add beta betaOffset
  updateStep maxIterations
             epsilon
             2 -- we've done one iteration so far, so the next one will be the second
             newBeta
             xOffset
             nrResults.sumLogLikelihood
             (updateTTEData newBeta xOffset tteData)

updateStep :: Int -> Double -> Int -> VS.Vector Double -> VS.Vector Double -> Double -> TTEData -> Either T.Text (VS.Vector Double)
updateStep maxIterations epsilon iteration beta xOffset logLikelihood tteData
  | iteration > maxIterations =
    Left "Did not converge in the alloted number of iterations"
  | otherwise =
    case calcBetaOffset tteData of
      Left e -> Left e
      Right (betaOffset, nrResults) ->
        let newBeta = L.add beta betaOffset
        in     -- Case: we've achieved convergence
            if | checkNRConvergence epsilon logLikelihood nrResults ->
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
