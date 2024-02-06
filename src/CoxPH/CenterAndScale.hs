{-# OPTIONS_GHC -Wno-name-shadowing #-}

module CoxPH.CenterAndScale (
    centerAndScaleCovs,
) where

import CoxPH.Data
import Data.Either (fromRight, isLeft)
import Data.Text qualified as T
import Data.Vector qualified as V
import Data.Vector.Storable qualified as VS
import Data.Vector.Unboxed qualified as VU

-- | Conditionally center and scale a design matrix
centerAndScaleCovs ::
    -- | The covariate vectors comprising the design matrix. The covariate vectors are required to all be the same length.
    V.Vector (VS.Vector Double) ->
    -- | The subject weights. It is required that the length of this vector is equal to the length of each of the vectors in the design matrix.
    VS.Vector Double ->
    -- | Indicators for whether each of the covariates should be centered and scaled. It is required that the length of this vector is equal to the number of vectors in the design matrix.
    V.Vector ScaleCovariateIndicator ->
    Either T.Text (V.Vector (VS.Vector Double, Double))
-- The error handling could be improved to (i) list all of the errors and (ii)
-- to list the indices of the failing vectors
centerAndScaleCovs xDesignMatrix weights scaleIndicators =
    let nCovs = V.length xDesignMatrix
        nScaleIndicators = V.length scaleIndicators
        sumWeights = VS.sum weights
     in if nCovs /= nScaleIndicators
            then Left "The length of first input is not the same as the length of the third input"
            else do
                weightedMeans <- calcWeightedMeans sumWeights xDesignMatrix weights scaleIndicators
                let covTuples =
                        V.zip5
                            xDesignMatrix
                            (V.replicate nCovs weights)
                            (V.convert weightedMeans)
                            (V.replicate nCovs sumWeights)
                            scaleIndicators

                traverse conditionallyCenterAndScaleCov covTuples
  where
    conditionallyCenterAndScaleCov (x, _, _, _, ScaleCovariateNo) =
        Right (x, 1)
    conditionallyCenterAndScaleCov (x, weights, mean, sumWeights, ScaleCovariateYes) =
        centerAndScaleCov x weights mean sumWeights

centerAndScaleCov ::
    VS.Vector Double ->
    VS.Vector Double ->
    Double ->
    Double ->
    Either T.Text (VS.Vector Double, Double)
centerAndScaleCov x weights mean sumWeights
    | VS.length x /= VS.length weights =
        Left $
            T.concat
                [ "The length of the first input (the covariate values) "
                , "is not the same as the length of the second input "
                , "(the weights)"
                ]
    | otherwise =
        let xCentered = VS.map (subtract mean) x
            covariatePairs = VU.zip (VU.convert xCentered) (VU.convert weights)
            weightedAbsCovSum = VU.foldl calcAbsProd 0 covariatePairs
         in if weightedAbsCovSum == 0
                then Left "Constant column"
                else
                    let scaleVal = sumWeights / weightedAbsCovSum
                     in Right (VS.map (* scaleVal) xCentered, scaleVal)
  where
    calcAbsProd :: Double -> (Double, Double) -> Double
    calcAbsProd sumVal (xVal, weightVal) = sumVal + abs (xVal * weightVal)

-- The error handling could be improved to (i) list all of the errors and (ii)
-- to list the indices of the failing vectors
calcWeightedMeans ::
    Double ->
    V.Vector (VS.Vector Double) ->
    VS.Vector Double ->
    V.Vector ScaleCovariateIndicator ->
    Either T.Text (VS.Vector Double)
calcWeightedMeans sumWeights xDesignMatrix weights scaleIndicators = do
    let covariatePairs = V.zip xDesignMatrix scaleIndicators
        eitherMeans =
            V.map
                (conditionallyCalcWeightedMean weights)
                covariatePairs
        errorCovs = V.filter isLeft eitherMeans
    if V.length errorCovs >= 1
        then -- The conditional statement makes the `Right` case unreachable
        case V.head errorCovs of
            Left x -> Left x
            Right x -> Right (VS.singleton x)
        else Right (VS.convert (V.map (fromRight 0) eitherMeans))
  where
    conditionallyCalcWeightedMean ::
        VS.Vector Double ->
        (VS.Vector Double, ScaleCovariateIndicator) ->
        Either T.Text Double
    conditionallyCalcWeightedMean _ (_, ScaleCovariateNo) = Right 0
    conditionallyCalcWeightedMean weights (x, ScaleCovariateYes) =
        calcWeightedMean x weights sumWeights

calcWeightedMean ::
    VS.Vector Double ->
    VS.Vector Double ->
    Double ->
    Either T.Text Double
calcWeightedMean x weights sumWeights
    | VS.length x == 0 = Left "Can't take the mean of a length-0 vector"
    | sumWeights == 0 = Left "Can't take the mean when the weights sum to 0"
    | otherwise =
        let sumCovs = VU.foldl op 0 (VU.zip (VU.convert x) (VU.convert weights))
         in Right (sumCovs / sumWeights)
  where
    op :: Double -> (Double, Double) -> Double
    op sumCovs (covVal, weightVal) = sumCovs + covVal * weightVal

-- -- Possibly don't need ---------------------------------------------------------

-- calcWeightedMeans2 :: V.Vector (Vector Double)
--                   -> Vector Double
--                   -> V.Vector ScaleCovariateIndicator
--                   -> Either T.Text (Vector Double)
-- calcWeightedMeans2 xDesignMatrix weights scaleIndicators = do
--   let covariatePairs = V.zip xDesignMatrix scaleIndicators
--   let vecEitherMeans = V.map (conditionallyCalcMean weights) covariatePairs
--   if V.any isRight vecEitherMeans
--   then Left "There was an error"  -- TODO: need to make this better
--   else Right (VS.convert (V.map (fromRight 0) vecEitherMeans))
--   where
--     conditionallyCalcMean :: Vector Double
--                           -> (Vector Double, ScaleCovariateIndicator)
--                           -> Either T.Text Double
--     conditionallyCalcMean _ (_, ScaleCovariateNo) = Right 0
--     conditionallyCalcMean weights (x, ScaleCovariateYes) = calcWeightedMean2 x weights

-- calcWeightedMean2 :: Vector Double
--                  -> Vector Double
--                  -> Either T.Text Double
-- calcWeightedMean2 x weights =
--   if VS.length x == 0
--   then Left "Can't take the mean of a length-0 vector"
--   else let (sumCovs, sumWeights) = VS.foldl op (0, 0) (VS.zip x weights)
--        in  if sumWeights == 0
--            then Left "Can't take the mean when the weights sum to 0"
--            else Right (sumCovs / sumWeights)
--   where
--     op :: (Double, Double) -> (Double, Double) -> (Double, Double)
--     op (sumCovs, sumWeights) (covVal, weightVal) =
--       (sumCovs + covVal * weightVal, sumWeights + weightVal)

-- calcUnsafeVecMean :: Vector Double -> Double
-- calcUnsafeVecMean x =
--   VS.sum x / fromIntegral (VS.length x)
