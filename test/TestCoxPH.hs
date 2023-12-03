{-# OPTIONS_GHC -Wno-name-shadowing #-}

module Main (main) where

import Control.Monad ( unless )
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import Data.Text (Text)
import Test.Tasty ( TestTree, defaultMain, testGroup )
import Test.Tasty.HUnit ( Assertion, testCase, (@?=), assertFailure )
import CoxPH

main :: IO ()
main = defaultMain unitTests

unitTests :: TestTree
unitTests = testGroup "Unit tests"
  [ testCase "Constant weights"
             (assertEqualCovs actualConstWeights expectedConstWeights)
  , testCase "Non-constant weights"
             (assertEqualCovs actualNonConstWeights expectedNonConstWeights)
  ]


-- Helper functions ------------------------------------------------------------

-- How much tolerance is permitted for floating point comparisons
eps :: Double
eps = 0.0000000001

-- Compare two objects of the same shape as created by `centerAndScaleCovs`
assertEqualCovs :: Either Text (V.Vector (VU.Vector Double, Double))
                -> Either Text (V.Vector (VU.Vector Double, Double))
                -> Assertion
assertEqualCovs (Left _) (Right _) = assertFailure "Should be a `Right`"
assertEqualCovs (Right _) (Left _) = assertFailure "Should be a `Left`"
assertEqualCovs (Left actual) (Left expected) = actual @?= expected
assertEqualCovs (Right actual) (Right expected)
  | length actual /= length expected = assertFailure "Top-level lengths don't match"
  | otherwise =
    if not (V.all (uncurry checkEqualScales) (V.zip actual expected))
    then assertFailure "Not all scales match"
    else unless (V.all (uncurry checkEqualCovs) (V.zip actual expected))
                (assertFailure "Not all covariates match")
  where
    checkEqualScales :: (a1, Double) -> (a2, Double) -> Bool
    checkEqualScales (_, v1) (_, v2) = approxEqual eps v1 v2
    checkEqualCovs :: (VU.Vector Double, a1) -> (VU.Vector Double, a2) -> Bool
    checkEqualCovs (x1, _) (x2, _)
      | VU.length x1 /= VU.length x2 = False
      | otherwise = VU.all (uncurry (approxEqual eps)) (VU.zip x1 x2)

-- Compare two numbers for equality up to the specified level of precision
approxEqual :: Double -> Double -> Double -> Bool
approxEqual eps x y = abs (x - y) <= eps


-- Create test data ------------------------------------------------------------

x1 :: VU.Vector Double
x1 = VU.fromList [1, 2, 3, 4, 5, 6]

x2 :: VU.Vector Double
x2 = VU.fromList [0, 1, 1, 0, 0, 2]

covs :: V.Vector (VU.Vector Double)
covs = V.fromList [x1, x2]

constWeights :: VU.Vector Double
constWeights = VU.fromList [1, 1, 1, 1, 1, 1]

nonConstWeights :: VU.Vector Double
nonConstWeights = VU.fromList [1, 2, 3, 4, 5, 6]

scaleIndicators :: V.Vector ScaleCovariateIndicator
scaleIndicators = V.fromList [ScaleCovariateYes, ScaleCovariateNo]


-- Create tests ----------------------------------------------------------------

actualConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
actualConstWeights = centerAndScaleCovs covs constWeights scaleIndicators

expectedConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
expectedConstWeights =
  let newX1 = VU.fromList
        [ -1.6666666666666665
        , -1
        , -0.3333333333333333
        ,  0.3333333333333333
        ,  1
        ,  1.6666666666666665
        ]
      expected = V.fromList
        [ (newX1, 0.6666666666666666)
        , (x2, 1)
        ]
  in  Right expected

actualNonConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
actualNonConstWeights = centerAndScaleCovs covs nonConstWeights scaleIndicators

expectedNonConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
expectedNonConstWeights =
  let newX1 = VU.fromList
        [ -2.6250000000000000
        , -1.8374999999999999
        , -1.0499999999999998
        , -0.2624999999999998
        ,  0.5250000000000002
        ,  1.3125000000000004
        ]
      expected = V.fromList
        [ (newX1, 0.7875000000000001 )
        , (x2, 1)
        ]
  in  Right expected
