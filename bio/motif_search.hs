import Data.List

bruteforceMedian = Nothing

dnaAlphabet = "ACTG"

-- Hamming distance
dH seq1 seq2 = sum $ zipWith dist seq1 seq2
                where dist a b | a==b      = 0
                               | otherwise = 1

total_dH l_mer seqs = Nothing

l_mers xs l = map (take l) $ take (n-l+1) $ tails xs
              where n = length xs