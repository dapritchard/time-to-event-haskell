module CoxPH (
  CoxPHConvergenceFailure(..),
  CoxPHMethod(..),
  Delta(..),
  LastInStrataIndicator(..),
  ScaleCovariateIndicator(..),
  coxph
  ) where

import qualified Data.Vector.Unboxed as VU

data Delta = ObservedEvent | Censored

-- type Vec2DimDouble = VU.Vector (VU.Vector Double)

data CoxPHConvergenceFailure = CoxPHConvergenceFailure

data CoxPHMethod = Breslow | Effron

data LastInStrataIndicator = LastInStrataNo | LastInStrataYes

data ScaleCovariateIndicator = ScaleCovariateNo | ScaleCovariateYes

-- type CoxPHResult = Either CoxPHConvergenceFailure (VU.Vector Double)
type CoxPHResult = Either Text (VU.Vector Double)

coxph :: VU.Vector Double                   -- length n, the per-subject minumum times of the event or censoring times
      -> VU.Vector Delta                    -- length n, indicators for whether patient observered an event or censoring
      -> VU.Vector (VU.Vector Double)        -- dimension n by p, X design matrix
      -> VU.Vector Double                   -- length n, X offset vector
      -> VU.Vector Double                   -- length n, subject weights
      -> VU.Vector LastInStrataIndicator    -- length n, indicators for whether a patient had the latest minumum event of censoring time within a strata
      -> VU.Vector Double                   -- length p, starting values for the beta coefficients
      -> VU.Vector ScaleCovariateIndicator  -- length p, indicators for whether a given covariate should be centered and scaled
      -> CoxPHMethod                       -- which of the Breslow or Effron methods to use in the event of tied event times
      -> Integer                           -- the maximumn number of Newton-Rhapson iterations to perform
      -> Double                            -- the value for which the absolute value of 1 minus ratio of the likelihood for two consecutive iterations must be below for convergence to be achieved
      -> Double                            -- TODO
      -> CoxPHResult
coxph _ -- time
      _ -- eventStatus
      _ -- xDesignMatrix
      _ -- xOffset
      _ -- weights
      _ -- strataIndicator
      _ -- beta
      _ -- scaleIndicator
      _ -- method
      _ -- maxIterations
      _ -- epsilon
      _ = do -- tolerance
  Right (VU.fromList [])

calcVecMean :: V.Vector Double -> Either CoxPHConvergenceFailure Double
calcVecMean x =
  if lengthX == 0
  then Left CoxPHConvergenceFailure
  else Right (V.sum x / fromIntegral lengthX)
  where
    lengthX = V.length x
