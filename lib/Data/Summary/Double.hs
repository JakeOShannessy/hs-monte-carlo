{-# LANGUAGE BangPatterns #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Summary.Double
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Summary statistics for @Double@s.
--

module Data.Summary.Double (
    -- * The @Summary@ data type
    Summary,
    summary,
    summaryH,
    Histogram (Histogram),
    Bin (Bin),
    getHist,
    sizeHist,
    update,

    -- * @Summary@ properties
    sampleSize,
    sampleMin,
    sampleMax,
    sampleMean,
    sampleSE,
    sampleVar,
    sampleSD,
    sampleCI,

    ) where

import Control.DeepSeq
import Data.List( foldl' )
import Data.Monoid
import Text.Printf

import Data.Summary.Utils


-- | A type for storing summary statistics for a data set including
-- sample size, min and max values, and first and second moments.
data Summary = S {-# UNPACK #-} !Int     -- sample size
                 {-# UNPACK #-} !Double  -- sample mean
                 {-# UNPACK #-} !Double  -- sum of squares
                 {-# UNPACK #-} !Double  -- sample min
                 {-# UNPACK #-} !Double  -- sample max
                 !(Maybe Histogram)              -- histogram of the sample data (unnormalised)
                 
getHist (S _ _ _ _ _ (Just hist)) = hist

data Histogram = Histogram
    !(Double, Double)    -- ^Upper and lower bounds
    !(Int, Int)  -- ^No. of out of bound values (upper, lower)
    ![Bin]   -- ^Bins
    deriving Show
    
sizeHist (Histogram _ _ bins) = foldl' (\acc (Bin _ _ n)->acc+n) 0 bins
    
data Bin = Bin
    !Double -- ^Lower limit inclusive
    !Double -- ^Upper limit exclusive
    !Int    -- ^Number of samples
    deriving Show
    
-- TODO: run the simulation for a certain number of samples to obtain reasonable histogram values
emptyHist :: Int -> (Double, Double) -> Histogram
emptyHist n (rangeLow, rangeHigh) = Histogram (minX, maxX) (0,0) bins
    where   -- TODO: Find a more thorough way of determining an appropriate range, removing outliers etc
        -- totalNo = fromIntegral (length sample) :: Double
        minX = rangeLow
        maxX = rangeHigh
        range = maxX - minX
        
        binSize = range/(fromIntegral n)   -- make sure all of this works well with floating point arithmetic
        -- bins = [minX,minX+binSize..maxX] -- it doesn't at the moment
        bins = mkBins [Bin (maxX-binSize) maxX 0]
            where
                mkBins ((Bin lo hi n):ps) | lo > minX = mkBins $  (Bin (lo-binSize) (hi-binSize) 0):(Bin lo hi n):ps
                                          | otherwise = ((Bin lo hi 0):ps)

instance Show Summary where
    show s@(S n mu _ l h hist) =
        printf "    sample size: %d" n
        ++ printf "\n            min: %g" l
        ++ printf "\n            max: %g" h
        ++ printf "\n           mean: %g" mu
        ++ printf "\n             SE: %g" (sampleSE s)
        ++ printf "\n         99%% CI: (%g, %g)" c1 c2
        ++ show hist
      where (c1,c2) = sampleCI 0.99 s

instance Monoid Summary where
    mempty = empty
    mappend = union

instance NFData Summary
instance NFData Histogram
instance NFData Bin

-- | Get a summary of a list of values.
summary :: [Double] -> Summary
summary = foldl' update empty

-- | Get a summary of a list of values with a histogram.
summaryH :: Int -> (Double, Double) -> [Double] -> Summary
summaryH n (min,max) = foldl' update $ emptyH n (min,max)

-- | Get an empty summary.
empty :: Summary
empty = S 0 0 0 (1/0) (-1/0) Nothing

-- | With a certain histogram.
emptyH :: Int -> (Double, Double) -> Summary
emptyH n (min,max) = S 0 0 0 (1/0) (-1/0) $ Just $ emptyHist n (min,max)

-- | Update the summary with a data point.
-- Running mean and variance computed as in Knuth, Vol 2, page 232,
-- 3rd edition, see http://www.johndcook.com/standard_deviation.html for
-- a description.
update :: Summary -> Double -> Summary
update (S n m s l h hist) x =
    let n'    = n+1
        delta = x - m
        m'    = m + delta / fromIntegral n'
        s'    = s + delta*(x - m')
        l'    = if x < l then x else l
        h'    = if x > h then x else h
        hist' = hist >>= (\hist -> return $ updateHist hist x)
    in S n' m' s' l' h' (deepseq hist' hist')
    
    
updateHist :: Histogram -> Double -> Histogram
updateHist (Histogram (lowRange, highRange) (tooLow, tooHigh) bins) value
    | value < lowRange = Histogram (lowRange, highRange) (tooLow+1, tooHigh) bins
    | value >= highRange = Histogram (lowRange, highRange) (tooLow, tooHigh+1) bins
    | otherwise = Histogram (lowRange, highRange) (tooLow, tooHigh) (updateHist' [] bins value)
    where
        inBin :: Bin -> Double -> Bool
        inBin (Bin lo hi n) val = val >= lo && val < hi
        updateHist' :: [Bin] -> [Bin] -> Double -> [Bin]
        -- updateHist' checked (b:bins) x = (Bin lo hi (n+1)):bins
        updateHist' checked (b:bins) x | inBin b x = checked ++ ((Bin lo hi (n+1)):bins)
                                       | otherwise = updateHist' (checked ++ [b]) bins x
            where
                (Bin lo hi n) = b
unionHist
    (Histogram (lowRangeA, highRangeA) (tooLowA, tooHighA) binsA)
    (Histogram (lowRangeB, highRangeB) (tooLowB, tooHighB) binsB)
    | checkSimilarity = Histogram
                            (lowRangeA, highRangeA)
                            (tooLowA + tooLowB, tooHighA + tooHighB)
                            (unionBins binsA binsB)
    | otherwise = error "histograms do not match"
    where
        checkSimilarity = lowRangeA == lowRangeB && highRangeA == highRangeB && length binsA == length binsB
                              
unionBins = zipWith unionBin
unionBin :: Bin -> Bin -> Bin
unionBin binA@(Bin loA hiA nA) binB@(Bin loB hiB nB)
    | checkSimilarity = Bin loA hiA (nA+nB)
    | otherwise       = error "bins do not match"
    where
        checkSimilarity = loA == loB && hiA == hiB

-- | Take the union of two summaries.
-- Use the updating rules from Chan et al. "Updating Formulae and a Pairwise
--   Algorithm for Computing Sample Variances," available at
-- ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf
union :: Summary -> Summary -> Summary
union (S na ma sa la ha hista) (S nb mb sb lb hb histb) =
    let delta = mb - ma
        (na', nb') = (fromIntegral na, fromIntegral nb)
        n  = na + nb
        n' = fromIntegral n
        weightedDelta = delta*nb'/n'
        m  | n == 0    = 0
           | otherwise = ma + weightedDelta
        s  | n == 0    = 0
           | otherwise = sa + sb + delta*na'*weightedDelta
        l  = min la lb
        h  = max ha hb
        hist = do
            a' <- hista
            b' <- histb
            return $ unionHist a' b'
    in S n m s l h hist

-- | Get the sample size.
sampleSize :: Summary -> Int
sampleSize (S n _ _ _ _ _) = n

-- | Get the sample mean.
sampleMean :: Summary -> Double
sampleMean (S _ m _ _ _ _) = m

-- | Get the sample variance.
sampleVar :: Summary -> Double
sampleVar (S n _ s _ _ _) = s / fromIntegral (n - 1)

-- | Get the sample standard deviation.
sampleSD :: Summary -> Double
sampleSD s = sqrt (sampleVar s)

-- | Get the sample standard error.
sampleSE :: Summary -> Double
sampleSE s = sqrt (sampleVar s / fromIntegral (sampleSize s))

-- | Get a Central Limit Theorem-based confidence interval for the mean
-- with the specified coverage level.  The level must be in the range @(0,1)@.
sampleCI :: Double -> Summary -> (Double,Double)
sampleCI level s = interval level (sampleMean s) (sampleSE s)

-- | Get the minimum of the sample.
sampleMin :: Summary -> Double
sampleMin (S _ _ _ l _ _) = l

-- | Get the maximum of the sample.
sampleMax :: Summary -> Double
sampleMax (S _ _ _ _ h _) = h
