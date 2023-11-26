module CoxPH (
  CoxPHConvergenceFailure(..),
  CoxPHMethod(..),
  Delta(..),
  LastInStrataIndicator(..),
  ScaleCovariateIndicator(..),
  coxph
  ) where

import qualified Data.Vector.Unboxed as V

data Delta = ObservedEvent | Censored

-- type Vec2DimDouble = V.Vector (V.Vector Double)

data CoxPHConvergenceFailure = CoxPHConvergenceFailure

data CoxPHMethod = Breslow | Effron

data LastInStrataIndicator = LastInStrataNo | LastInStrataYes

data ScaleCovariateIndicator = ScaleCovariateNo | ScaleCovariateYes

type CoxPHResult = Either CoxPHConvergenceFailure (V.Vector Double)

coxph :: V.Vector Double                   -- length n, the per-subject minumum times of the event or censoring times
      -> V.Vector Delta                    -- length n, indicators for whether patient observered an event or censoring
      -> V.Vector (V.Vector Double)        -- dimension n by p, X design matrix
      -> V.Vector Double                   -- length n, X offset vector
      -> V.Vector Double                   -- length n, subject weights
      -> V.Vector LastInStrataIndicator    -- length n, indicators for whether a patient had the latest minumum event of censoring time within a strata
      -> V.Vector Double                   -- length p, starting values for the beta coefficients
      -> V.Vector ScaleCovariateIndicator  -- length p, indicators for whether a given covariate should be centered and scaled
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
  Right (V.fromList [])

calcVecMean :: V.Vector Double -> Either CoxPHConvergenceFailure Double
calcVecMean x =
  if lengthX == 0
  then Left CoxPHConvergenceFailure
  else Right (V.sum x / fromIntegral lengthX)
  where
    lengthX = V.length x
