-- |

module TestCoxPH.TestCoxPHUpdateNewtonRaphson where

import CoxPH.CenterAndScale
import CoxPH.CoxPHUpdateNewtonRaphson
import CoxPH.Data
import Data.Vector qualified as V
import Data.Vector.Storable qualified as VS
import Data.Either (fromRight)
-- import Numeric.LinearAlgebra

-- Create test input data ------------------------------------------------------

-- The following are the data inputs resulting from the following call
-- immediately before being input to `Ccoxfit6_iter`
--
-- test1 <- list(time=c(4,3,1,1,2,2,3),
--                status=c(1,1,1,0,1,1,0),
--                x=c(0,2,1,1,1,0,0),
--                sex=c(0,0,0,0,1,1,1))
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
times = VS.fromList [1, 1, 2, 2, 3, 3, 4]
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
beta :: VS.Vector Double
beta = VS.fromList [0, 0]
scaleIndicators :: V.Vector ScaleCovariateIndicator
scaleIndicators = V.fromList [ScaleCovariateYes, ScaleCovariateNo]
tiesMethod :: CoxPHMethod
tiesMethod = Breslow
maxIterations :: Int
maxIterations = 20
epsilon :: Double
epsilon = 0.000000000001818989

centeredAndScaledCovsResults :: V.Vector (VS.Vector Double, Double)
centeredAndScaledCovsResults =
  fromRight (V.singleton (VS.singleton 0, 0))
            (centerAndScaleCovs xDesignDataFrame
                                testWeights
                                scaleIndicators)
