{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module CoxPH (
  CoxPHConvergenceFailure(..),
  CoxPHMethod(..),
  Delta(..),
  LastInStrataIndicator(..),
  ScaleCovariateIndicator(..),
  centerAndScaleCovs,
  coxph
  ) where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import Data.Either (fromRight, isLeft)
import qualified Data.Text as T

data Delta = ObservedEvent | Censored

-- type Vec2DimDouble = VU.Vector (VU.Vector Double)

data CoxPHConvergenceFailure = CoxPHConvergenceFailure

data CoxPHMethod = Breslow | Effron

data LastInStrataIndicator = LastInStrataNo | LastInStrataYes

data ScaleCovariateIndicator = ScaleCovariateNo | ScaleCovariateYes

-- type CoxPHResult = Either CoxPHConvergenceFailure (VU.Vector Double)
type CoxPHResult = Either T.Text (VU.Vector Double)

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

-- The error handling could be improved to (i) list all of the errors and (ii)
-- to list the indices of the failing vectors
centerAndScaleCovs :: V.Vector (VU.Vector Double)
                   -> VU.Vector Double
                   -> V.Vector ScaleCovariateIndicator
                   -> Either T.Text (V.Vector (VU.Vector Double, Double))
centerAndScaleCovs xDesignMatrix weights scaleIndicators =
  let nCovs = V.length xDesignMatrix
      nScaleIndicators = V.length scaleIndicators
      sumWeights = VU.sum weights
  in if nCovs /= nScaleIndicators
     then Left "The length of first input is not the same as the length of the third input"
     else do
       weightedMeans <- calcWeightedMeans sumWeights xDesignMatrix weights scaleIndicators
       let covTuples = V.zip5 xDesignMatrix
                              (V.replicate nCovs weights)
                              (V.convert weightedMeans)
                              (V.replicate nCovs sumWeights)
                              scaleIndicators
           eitherResults = V.map conditionallyCenterAndScaleCov covTuples
       if V.any isLeft eitherResults
       then Left "Not all covariates where able to be centered and scaled"
       else Right $ V.map (fromRight (VU.singleton 0, 0)) eitherResults
  where
    conditionallyCenterAndScaleCov (x, _, _, _, ScaleCovariateNo) =
      Right (x, 1)
    conditionallyCenterAndScaleCov (x, weights, mean, sumWeights, ScaleCovariateYes) =
      centerAndScaleCov x weights mean sumWeights

centerAndScaleCov :: VU.Vector Double
                  -> VU.Vector Double
                  -> Double
                  -> Double
                  -> Either T.Text (VU.Vector Double, Double)
centerAndScaleCov x weights mean sumWeights
  | VU.length x /= VU.length weights =
      Left $ T.concat [ "The length of the first input (the covariate values) "
                      , "is not the same as the length of the second input "
                      , "(the weights)"
                      ]
  | otherwise =
      let xCentered = VU.map (subtract mean) x
          covariatePairs = VU.zip xCentered weights
          weightedAbsCovSum = VU.foldl calcAbsProd 0 covariatePairs
      in  if weightedAbsCovSum == 0
          then Left "Constant column"
          else let scaleVal = sumWeights / weightedAbsCovSum
               in  Right (VU.map (* scaleVal) xCentered, scaleVal)
  where
    calcAbsProd :: Double -> (Double, Double) -> Double
    calcAbsProd sumVal (xVal, weightVal) = sumVal + abs (xVal * weightVal)

-- The error handling could be improved to (i) list all of the errors and (ii)
-- to list the indices of the failing vectors
calcWeightedMeans :: Double
                  -> V.Vector (VU.Vector Double)
                  -> VU.Vector Double
                  -> V.Vector ScaleCovariateIndicator
                  -> Either T.Text (VU.Vector Double)
calcWeightedMeans sumWeights xDesignMatrix weights scaleIndicators = do
  let covariatePairs = V.zip xDesignMatrix scaleIndicators
      eitherMeans = V.map (conditionallyCalcWeightedMean weights)
                             covariatePairs
      errorCovs = V.filter isLeft eitherMeans
  if V.length errorCovs >= 1
  then
    -- The conditional statement makes the `Right` case unreachable
    case V.head errorCovs of
      Left x -> Left x
      Right x -> Right (VU.singleton x)
  else Right (VU.convert (V.map (fromRight 0) eitherMeans))
  where
    conditionallyCalcWeightedMean :: VU.Vector Double
                                  -> (VU.Vector Double, ScaleCovariateIndicator)
                                  -> Either T.Text Double
    conditionallyCalcWeightedMean _ (_, ScaleCovariateNo) = Right 0
    conditionallyCalcWeightedMean weights (x, ScaleCovariateYes) =
      calcWeightedMean x weights sumWeights

calcWeightedMean :: VU.Vector Double
                 -> VU.Vector Double
                 -> Double
                 -> Either T.Text Double
calcWeightedMean x weights sumWeights
  | VU.length x == 0 = Left "Can't take the mean of a length-0 vector"
  | sumWeights == 0 = Left "Can't take the mean when the weights sum to 0"
  | otherwise =
    let sumCovs = VU.foldl op 0 (VU.zip x weights)
    in  Right (sumCovs / sumWeights)
  where
    op :: Double -> (Double, Double) -> Double
    op sumCovs (covVal, weightVal) = sumCovs + covVal * weightVal


-- -- Possibly don't need ---------------------------------------------------------

-- calcWeightedMeans2 :: V.Vector (VU.Vector Double)
--                   -> VU.Vector Double
--                   -> V.Vector ScaleCovariateIndicator
--                   -> Either T.Text (VU.Vector Double)
-- calcWeightedMeans2 xDesignMatrix weights scaleIndicators = do
--   let covariatePairs = V.zip xDesignMatrix scaleIndicators
--   let vecEitherMeans = V.map (conditionallyCalcMean weights) covariatePairs
--   if V.any isRight vecEitherMeans
--   then Left "There was an error"  -- TODO: need to make this better
--   else Right (VU.convert (V.map (fromRight 0) vecEitherMeans))
--   where
--     conditionallyCalcMean :: VU.Vector Double
--                           -> (VU.Vector Double, ScaleCovariateIndicator)
--                           -> Either T.Text Double
--     conditionallyCalcMean _ (_, ScaleCovariateNo) = Right 0
--     conditionallyCalcMean weights (x, ScaleCovariateYes) = calcWeightedMean2 x weights

-- calcWeightedMean2 :: VU.Vector Double
--                  -> VU.Vector Double
--                  -> Either T.Text Double
-- calcWeightedMean2 x weights =
--   if VU.length x == 0
--   then Left "Can't take the mean of a length-0 vector"
--   else let (sumCovs, sumWeights) = VU.foldl op (0, 0) (VU.zip x weights)
--        in  if sumWeights == 0
--            then Left "Can't take the mean when the weights sum to 0"
--            else Right (sumCovs / sumWeights)
--   where
--     op :: (Double, Double) -> (Double, Double) -> (Double, Double)
--     op (sumCovs, sumWeights) (covVal, weightVal) =
--       (sumCovs + covVal * weightVal, sumWeights + weightVal)

-- calcUnsafeVecMean :: VU.Vector Double -> Double
-- calcUnsafeVecMean x =
--   VU.sum x / fromIntegral (VU.length x)
