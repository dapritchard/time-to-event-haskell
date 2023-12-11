{-# OPTIONS_GHC -Wno-name-shadowing #-}

module CoxPH.CenterAndScale (
  centerAndScaleCovs
  ) where

import CoxPH.Data
import Data.Vector.Unboxed qualified as VU
import Data.Vector qualified as V
import Data.Text qualified as T
import Data.Either (fromRight, isLeft)

{- | Conditionally center and scale a design matrix
-}
centerAndScaleCovs :: V.Vector (VU.Vector Double)
                   -- ^ The covariate vectors comprising the design matrix. The covariate vectors are required to all be the same length.
                   -> VU.Vector Double
                   -- ^ The subject weights. It is required that the length of this vector is equal to the length of each of the vectors in the design matrix.
                   -> V.Vector ScaleCovariateIndicator
                   -- ^ Indicators for whether each of the covariates should be centered and scaled. It is required that the length of this vector is equal to the number of vectors in the design matrix.
                   -> Either T.Text (V.Vector (VU.Vector Double, Double))
-- The error handling could be improved to (i) list all of the errors and (ii)
-- to list the indices of the failing vectors
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
