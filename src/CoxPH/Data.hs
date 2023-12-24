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
  deriving Show

data CoxPHConvergenceFailure = CoxPHConvergenceFailure
  deriving Show

data CoxPHMethod = Breslow | Efron
  deriving Show

data LastInStrataIndicator = LastInStrataNo | LastInStrataYes
  deriving Show

data ScaleCovariateIndicator = ScaleCovariateNo | ScaleCovariateYes
  deriving Show
