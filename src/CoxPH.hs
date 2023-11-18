module CoxPH (
  Delta(..),
  coxph
  ) where

import qualified Data.Vector.Unboxed as V

data Delta = ObservedEvent | Censored

-- type Vec2DimDouble = V.Vector (V.Vector Double)

newtype CoxPHResult = CoxPHResult (V.Vector Double)

coxph :: V.Vector Double -> V.Vector Delta -> V.Vector (V.Vector Double) -> CoxPHResult
coxph _ _ _ = CoxPHResult (V.fromList [])
