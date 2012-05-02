module MotifSearch where

import Control.Monad
import Data.List
import TreeUtils
import Data.Tree

-- Find motifs of length k in list of DNA sequences using "brute force" algorithm.
-- Searching for minimum total distance between sequences and all 4^k k-mers
bruteforceMedian seqs k = map snd $ filter f $ zip totalDistances words
                          where totalDistances = map (total_dH seqs) words
                                words =  replicateM k dnaAlphabet
                                f (n, _) = n == minimum totalDistances
dnaAlphabet = "ACTG"

-- Hamming distance
dH seq1 seq2 = sum $ zipWith dist seq1 seq2
                where dist a b | a==b      = 0
                               | otherwise = 1

-- Total Hamming Distance
total_dH seqs k_mer = sum $ map min_dH seqs
                      where min_dH seq = minimum $ map (dH k_mer) $ k_mers (length k_mer) seq

-- All k-mers from sequence
-- For exapmle: 
-- k_mers 3 "AACCTT" 
-- ["AA","AC","CC","CT","TT"]
k_mers k seq = map (take k) $ take (n-k+1) $ tails seq
              where n = length seq

-- Branch And Bound median string search

nextLevel k prefix | length prefix < k      = (prefix, map (appendChar prefix) dnaAlphabet)
                   | otherwise              = (prefix, [])
                     where appendChar str c = str ++ [c]

sampleSeqs = ["TGACCGTGCAAAAAA", "TAGAAGAAAAAATGG", "AAAAAAATCATGACT"]

data TreeItem = Item {k_mer :: String,
                      distance :: Int}
                deriving (Show, Eq)

instance Ord TreeItem where
      compare i1 i2 = compare (distance i1) (distance i2)

searchTree seqs k = mapTree f $ unfoldTree (nextLevel k) ""
                    where f x = Item x (total_dH seqs x)


step alts = minimum $ map rootLabel $ alts

availableVertices tree k | null nextToLastVertices = []
                         | otherwise = subForest $ head $ nextToLastVertices
                           where nextToLastVertices | null vertices  = [] 
                                                    | otherwise =  dropWhile (null . subForest) $ head vertices 
                                                      where vertices = drop (k-1) (levels' tree)

improve k sTree seq  | null nextVertices = seq
                     | otherwise  =  improve k (prune (>= minItem) sTree) (k_mer minItem)
                                    where minItem = step nextVertices
                                          nextVertices = availableVertices sTree k 


bbMedianString seqs k = improve k sTree ""
                        where sTree = searchTree seqs k

