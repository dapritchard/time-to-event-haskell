module CoxPH.OrderByStrata where

import Data.Vector qualified as V
import Data.Vector.Algorithms.Tim (sortBy)
import Data.Vector.Storable qualified as VS

orderByStrata :: VS.Vector Int -> VS.Vector Double -> VS.Vector Int
orderByStrata strata times =
    let orderedTuples = orderByStrataTuples strata times
     in VS.convert (V.map extractIndex orderedTuples)
  where
    extractIndex :: (Int, Int, Double) -> Int
    extractIndex (index, _, _) = index

orderByStrataTuples ::
    VS.Vector Int ->
    VS.Vector Double ->
    V.Vector (Int, Int, Double)
orderByStrataTuples strata times = V.create $ do
    let n = VS.length strata
        entriesOrig =
            V.zip3
                (V.iterateN n (+ 1) 0)
                (V.convert strata)
                (V.convert times)
    entriesMutable <- V.thaw entriesOrig
    sortBy compareEntries entriesMutable
    return entriesMutable
  where
    compareEntries :: (Int, Int, Double) -> (Int, Int, Double) -> Ordering
    compareEntries (_, stratum1, time1) (_, stratum2, time2)
        | stratum1 < stratum2 = LT
        | stratum1 == stratum2 = compare time1 time2
        | otherwise = GT
