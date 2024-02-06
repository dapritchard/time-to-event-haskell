module CoxPH (
    CoxPHConvergenceFailure (..),
    CoxPHMethod (..),
    Delta (..),
    LastInStrataIndicator (..),
    ScaleCovariateIndicator (..),
    coxph,
) where

import CoxPH.CenterAndScale (centerAndScaleCovs)
import CoxPH.CoxPHUpdateNewtonRaphson (NRResults (..), TTEData (..), coxPHUpdateNewtonRaphson)
import CoxPH.Data
import Data.Text qualified as T
import Data.Vector qualified as V
import Data.Vector.Storable qualified as VS
import Numeric.LinearAlgebra qualified as L

-- | TODO:
coxph ::
    -- | length n, the per-subject minumum times of the event or censoring times
    VS.Vector Double ->
    -- | length n, indicators for whether patient observered an event or censoring
    V.Vector Delta ->
    -- | dimension n by p, X design matrix
    V.Vector (VS.Vector Double) ->
    -- | length n, X offset vector
    VS.Vector Double ->
    -- | length n, subject weights
    VS.Vector Double ->
    -- | length n, strata
    VS.Vector Int ->
    -- | length p, starting values for the beta coefficients
    VS.Vector Double ->
    -- | length p, indicators for whether a given covariate should be centered and scaled
    V.Vector ScaleCovariateIndicator ->
    -- | which of the Breslow or Effron methods to use in the event of tied event times
    CoxPHMethod ->
    -- | the maximum number of Newton-Rhapson iterations to perform
    Int ->
    -- | the value for which the absolute value of 1 minus ratio of the likelihood for two consecutive iterations must be below for convergence to be achieved
    -- -> Double                            -- TODO
    Double ->
    Either T.Text (VS.Vector Double)
coxph
    times
    eventStatuses
    xDesignDataFrame
    xOffset
    weights
    strata
    beta
    scaleIndicators
    tiesMethod
    maxIterations
    epsilon =
        -- tolerance
        do
            centeredAndScaledCovsResults <-
                centerAndScaleCovs
                    xDesignDataFrame
                    weights
                    scaleIndicators
            let (centeredAndScaledCovs, scales) = V.unzip centeredAndScaledCovsResults
                -- scales = V.map snd centeredAndScaledCovsResults
                xDesignMatrix = L.fromColumns (V.toList centeredAndScaledCovs)
                tteData =
                    TTEData
                        { time = times
                        , eventStatus = eventStatuses
                        , xDesignMatrix = xDesignMatrix
                        , stratum = strata
                        , weights = weights
                        , xProdBeta = L.add xOffset (xDesignMatrix L.#> beta)
                        , tiesMethod = tiesMethod
                        }
            (betaOffset, nrResults) <- calcBetaOffset tteData
            let newBeta = L.add beta betaOffset
            finalBeta <-
                updateStep
                    maxIterations
                    epsilon
                    2 -- since we've manually done one iteration so far
                    newBeta
                    xOffset
                    nrResults.sumLogLikelihood
                    (updateTTEData newBeta xOffset tteData)
            Right (VS.convert scales / finalBeta)

updateStep :: Int -> Double -> Int -> VS.Vector Double -> VS.Vector Double -> Double -> TTEData -> Either T.Text (VS.Vector Double)
updateStep maxIterations epsilon iteration beta xOffset logLikelihood tteData
    | iteration > maxIterations =
        Left "Did not converge in the alloted number of iterations"
    | otherwise =
        case calcBetaOffset tteData of
            Left e -> Left e
            Right (betaOffset, nrResults) ->
                let newBeta = L.add beta betaOffset
                 in if checkNRConvergence epsilon logLikelihood nrResults
                        then -- Case: we've achieved convergence;
                            Right newBeta
                        else -- Case: we haven't achieved convergence yet, but the log
                             -- likelihood is moving in the right direction

                            updateStep
                                maxIterations
                                epsilon
                                (iteration + 1)
                                newBeta
                                xOffset
                                nrResults.sumLogLikelihood
                                (updateTTEData newBeta xOffset tteData)

updateTTEData :: VS.Vector Double -> VS.Vector Double -> TTEData -> TTEData
updateTTEData beta xOffset tteData =
    let newXProdBeta =
            L.add
                xOffset
                (tteData.xDesignMatrix L.#> beta)
     in tteData{xProdBeta = newXProdBeta}

-- calculates I^{-1}(\hat{\beta}^{(n)})\, U(\hat{\beta}^{(n)}), i.e. the term
-- that we add to \hat{\beta}^{(n)} to update it
calcBetaOffset :: TTEData -> Either T.Text (VS.Vector Double, NRResults)
calcBetaOffset tteData = do
    let nrResults = coxPHUpdateNewtonRaphson tteData
        inverseInformationMatrix = L.inv nrResults.informationMatrix -- FIXME: this can fail. Need to catch failure as appropriate
        betaOffset = inverseInformationMatrix L.#> nrResults.score
    Right (betaOffset, nrResults)

checkNRConvergence :: Double -> Double -> NRResults -> Bool
checkNRConvergence epsilon logLikelihood nrResults =
    abs (1 - (logLikelihood / nrResults.sumLogLikelihood)) < epsilon
