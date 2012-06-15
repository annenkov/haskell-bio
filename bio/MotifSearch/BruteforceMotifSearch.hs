module BruteforceMotifSearch where

import Control.Monad
import Data.List

-- Find motifs of length k in list of DNA sequences using "brute force" algorithm.
-- Searching for minimum total distance between sequences and all 4^k k-mers
bruteforceMedianString seqs k = map snd $ filter minDist $ zip totalDistances words
                          where totalDistances = map (total_dH seqs) words
                                words =  replicateM k dnaAlphabet
                                minDist (n, _) = n == minimum totalDistances

dnaAlphabet = "ACGT"

-- Hamming distance
dH seq1 seq2 = sum $ zipWith diff seq1 seq2
               where diff a b | a == b    = 0
                              | otherwise = 1

-- Total Hamming Distance
total_dH seqs k_mer = sum $ map (min_dH k_mer) seqs
                      where min_dH k_mer seq = minimum $ map (dH k_mer) $ k_mers (length k_mer) seq

-- All k-mers from sequence
-- For exapmle: 
-- k_mers 3 "AACCTT" 
-- ["AA","AC","CC","CT","TT"]
k_mers k seq = map (take k) $ take (n-k+1) $ tails seq
               where n = length seq
