module CoxPH (
  CoxPHConvergenceFailure(..),
  Delta(..),
  coxph
  ) where

import qualified Data.Vector.Unboxed as V

data Delta = ObservedEvent | Censored

-- type Vec2DimDouble = V.Vector (V.Vector Double)

data CoxPHConvergenceFailure = CoxPHConvergenceFailure

type CoxPHResult = Either CoxPHConvergenceFailure (V.Vector Double)

coxph :: V.Vector Double -> V.Vector Delta -> V.Vector (V.Vector Double) -> CoxPHResult
coxph _ _ _ = Right (V.fromList [])
