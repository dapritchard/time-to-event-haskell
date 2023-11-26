{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module CoxPH (
  CoxPHConvergenceFailure(..),
  CoxPHMethod(..),
  Delta(..),
  LastInStrataIndicator(..),
  ScaleCovariateIndicator(..),
  coxph
  ) where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
-- import Data.Maybe (fromMaybe)
import Data.Either (fromRight, isRight)

import Data.Text (Text)

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
      -> V.Vector (VU.Vector Double)        -- dimension n by p, X design matrix
      -> VU.Vector Double                   -- length n, X offset vector
      -> VU.Vector Double                   -- length n, subject weights
      -> VU.Vector LastInStrataIndicator    -- length n, indicators for whether a patient had the latest minumum event of censoring time within a strata
      -> VU.Vector Double                   -- length p, starting values for the beta coefficients
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
      _ -- strataIndicators
      _ -- beta
      _ -- scaleIndicators
      _ -- method
      _ -- maxIterations
      _ -- epsilon
      _ = do -- tolerance
  Right (VU.fromList [])

-- centerAndScale :: VU.Vector (VU.Vector Double)
--                -> VU.Vector ScaleCovariateIndicator
--                -> Either Text VU.Vector (VU.Vector Double)
-- centerAndScale xDesignMatrix scaleIndicators = do
--   let vecMaybeMeans = VU.map calcWeightedMean xDesignMatrix
--   vecMeans <- if VU.any isNothing vecMeans
--               then Left "Had a length-0 vector"
--               else
--   where
--     isNothing :: Maybe -> Bool
--     isNothing Nothing -> True
--     isNothing _ -> False

calcWeightedMeans :: Double
                  -> V.Vector (VU.Vector Double)
                  -> VU.Vector Double
                  -> V.Vector ScaleCovariateIndicator
                  -> Either Text (VU.Vector Double)
calcWeightedMeans sumWeights xDesignMatrix weights scaleIndicators = do
  let covariatePairs = V.zip xDesignMatrix scaleIndicators
  let vecEitherMeans = V.map (conditionallyCalcMean weights) covariatePairs
  if V.any isRight vecEitherMeans
  then Left "There was an error"  -- TODO: need to make this better
  else Right (VU.convert (V.map (fromRight 0) vecEitherMeans))
  where
    conditionallyCalcMean :: VU.Vector Double
                          -> (VU.Vector Double, ScaleCovariateIndicator)
                          -> Either Text Double
    conditionallyCalcMean _ (_, ScaleCovariateNo) = Right 0
    conditionallyCalcMean weights (x, ScaleCovariateYes) = calcWeightedMean x weights sumWeights

calcWeightedMean :: VU.Vector Double
                 -> VU.Vector Double
                 -> Double
                 -> Either Text Double
calcWeightedMean x weights sumWeights =
  if VU.length x == 0
  then Left "Can't take the mean of a length-0 vector"
  else if sumWeights == 0
       then Left "Can't take the mean when the weights sum to 0"
       else let sumCovs = VU.foldl op 0 (VU.zip x weights)
            in Right (sumCovs / sumWeights)
  where
    op :: Double -> (Double, Double) -> Double
    op sumCovs (covVal, weightVal) = sumCovs + covVal * weightVal


-- Possibly don't need ---------------------------------------------------------

calcWeightedMeans2 :: V.Vector (VU.Vector Double)
                  -> VU.Vector Double
                  -> V.Vector ScaleCovariateIndicator
                  -> Either Text (VU.Vector Double)
calcWeightedMeans2 xDesignMatrix weights scaleIndicators = do
  let covariatePairs = V.zip xDesignMatrix scaleIndicators
  let vecEitherMeans = V.map (conditionallyCalcMean weights) covariatePairs
  if V.any isRight vecEitherMeans
  then Left "There was an error"  -- TODO: need to make this better
  else Right (VU.convert (V.map (fromRight 0) vecEitherMeans))
  where
    conditionallyCalcMean :: VU.Vector Double
                          -> (VU.Vector Double, ScaleCovariateIndicator)
                          -> Either Text Double
    conditionallyCalcMean _ (_, ScaleCovariateNo) = Right 0
    conditionallyCalcMean weights (x, ScaleCovariateYes) = calcWeightedMean2 x weights

calcWeightedMean2 :: VU.Vector Double
                 -> VU.Vector Double
                 -> Either Text Double
calcWeightedMean2 x weights =
  if VU.length x == 0
  then Left "Can't take the mean of a length-0 vector"
  else let (sumCovs, sumWeights) = VU.foldl op (0, 0) (VU.zip x weights)
       in  if sumWeights == 0
           then Left "Can't take the mean when the weights sum to 0"
           else Right (sumCovs / sumWeights)
  where
    op :: (Double, Double) -> (Double, Double) -> (Double, Double)
    op (sumCovs, sumWeights) (covVal, weightVal) =
      (sumCovs + covVal * weightVal, sumWeights + weightVal)

calcUnsafeVecMean :: VU.Vector Double -> Double
calcUnsafeVecMean x =
  VU.sum x / fromIntegral (VU.length x)
