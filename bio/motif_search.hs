import Data.List
import Control.Monad

bruteforceMedian seqs k = filter f $ zip (map (total_dH seqs) $ words) words
                          where totalDistances = (map (total_dH seqs) $ words)                                              
                                words =  replicateM 2 dnaAlphabet
                                f (n, k_mer) = n == minimum totalDistances
dnaAlphabet = "ACTG"

-- Hamming distance
dH seq1 seq2 = sum $ zipWith dist seq1 seq2
                where dist a b | a==b      = 0
                               | otherwise = 1

-- Total Hamming Distance
total_dH seqs k_mer = sum $ map min_dH seqs
                      where min_dH seq = minimum $ map (dH k_mer) $ k_mers (length k_mer) seq

-- all k-mers from sequence
k_mers k seq = map (take k) $ take (n-k+1) $ tails seq
              where n = length seq
