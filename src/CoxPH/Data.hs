-- |

module CoxPH.Data (
    CoxPHConvergenceFailure(..)
  , CoxPHMethod(..)
  , Delta(..)
  , LastInStrataIndicator(..)
  , ScaleCovariateIndicator(..)
  )
where

data Delta = ObservedEvent | Censored

data CoxPHConvergenceFailure = CoxPHConvergenceFailure

data CoxPHMethod = Breslow | Efron

data LastInStrataIndicator = LastInStrataNo | LastInStrataYes

data ScaleCovariateIndicator = ScaleCovariateNo | ScaleCovariateYes
