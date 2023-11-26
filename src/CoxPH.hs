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

coxph :: V.Vector Double             -- the per-subject minumums of the event or censoring times
      -> V.Vector Delta              -- indicators for whether patient observered an event or censoring
      -> V.Vector (V.Vector Double)  -- X design matrix
      -> V.Vector Double             --
      -> V.Vector Double
      -> V.Vector LastInStrataIndicator
      -> V.Vector Double
      -> V.Vector ScaleCovariateIndicator
      -> CoxPHMethod
      -> Integer
      -> Double
      -> Double
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
      _ = -- tolerance
  Right (V.fromList [])
