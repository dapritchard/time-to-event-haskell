module CoxPH (
  CoxPHConvergenceFailure(..),
  CoxPHMethod(..),
  Delta(..),
  LastInStrataIndicator(..),
  ScaleCovariateIndicator(..),
  coxph
  ) where

import CoxPH.CenterAndScale (centerAndScaleCovs)
import CoxPH.CoxPHUpdateNewtonRaphson ( TTEData(..) )
import CoxPH.Data
import Data.Vector.Storable qualified as VS
import Data.Vector qualified as V
import Data.Text qualified as T
import Numeric.LinearAlgebra qualified as L

-- type CoxPHResult = Either CoxPHConvergenceFailure (VU.Vector Double)
type CoxPHResult = Either T.Text (VS.Vector Double)

coxph :: VS.Vector Double                   -- length n, the per-subject minumum times of the event or censoring times
      -> V.Vector Delta                    -- length n, indicators for whether patient observered an event or censoring
      -> V.Vector (VS.Vector Double)     -- dimension n by p, X design matrix
      -> VS.Vector Double                   -- length n, X offset vector
      -> VS.Vector Double                   -- length n, subject weights
      -> VS.Vector Int                      -- length n, indicators for whether a patient had the latest minumum event of censoring time within a strata
      -> VS.Vector Double                   -- length p, starting values for the beta coefficients
      -> V.Vector ScaleCovariateIndicator  -- length p, indicators for whether a given covariate should be centered and scaled
      -> CoxPHMethod                        -- which of the Breslow or Effron methods to use in the event of tied event times
      -> Integer                            -- the maximumn number of Newton-Rhapson iterations to perform
      -> Double                             -- the value for which the absolute value of 1 minus ratio of the likelihood for two consecutive iterations must be below for convergence to be achieved
      -> Double                             -- TODO
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
      _ -- maxIterations
      _ -- epsilon
      _ = do -- tolerance
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
  Right (VS.fromList [])
