import Data.List
import Control.Monad
import Control.Monad.Random
import Test.HUnit

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


-- Random sequences for testing --

rndNucleotide :: Rand StdGen Char
rndNucleotide = do i <- getRandomR (0,3)
                   return $ dnaAlphabet!!i

rndDNA n = sequence $ replicate n rndNucleotide

-- generate m DNA sequences of length n
rndSeqs m n = sequence $ replicate m $ rndDNA n

implantMotif k_mer seq = do dna <- seq
                            i <- getRandomR (0,length dna - length k_mer - 1)
                            return $ take i dna ++ k_mer ++ drop i dna

makeTestSeqs m n k_mer =  sequence $ replicate m $ implantMotif k_mer $ rndDNA n


-- Tests --

test1 = TestCase (assertEqual "Finding 4-mer:" [k_mer] (bruteforceMedian testData 4))
        where k_mer = "AAAA"
              testData = evalRand (makeTestSeqs 5 20 k_mer) (mkStdGen 1)

-- run tests with "runTestTT tests"
tests = TestList [test1]