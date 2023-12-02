module Main (main) where

import Control.Monad ( unless )
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
-- import Test.HUnit.Approx ( assertApproxEqual )
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

eps :: Double
eps = 0.0000000001

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

approxEqual :: Double -> Double -> Double -> Bool
approxEqual eps x y = abs (x - y) <= eps

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

actualConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
actualConstWeights = centerAndScaleCovs covs constWeights scaleIndicators

expectedConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
expectedConstWeights =
  let newX1 = VU.fromList
        [ -1.66666666666667
        , -1
        , -0.333333333333333
        , 0.333333333333333
        , 1
        , 1.66666666666667
        ]
      expected = V.fromList
        [ (newX1, 0.666666666666667)
        , (x2, 1)
        ]
  in  Right expected

actualNonConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
actualNonConstWeights = centerAndScaleCovs covs nonConstWeights scaleIndicators

expectedNonConstWeights :: Either Text (V.Vector (VU.Vector Double, Double))
expectedNonConstWeights =
  let newX1 = VU.fromList
        [ -2.625, -1.837, -1.050, -0.262,  0.525,  1.313]
      expected = V.fromList
        [ (newX1, 0.7875)
        , (x2, 1)
        ]
  in  Right expected
