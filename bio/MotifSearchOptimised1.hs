module MotifSearchOptimised where

import Control.Monad
import Data.List
import TreeUtils
import Data.Tree
import System (getArgs)
import Bio.Sequence.Fasta
import Bio.Sequence.SeqData

-- Find motifs of length k in list of DNA sequences using "brute force" algorithm.
-- Searching for minimum total distance between sequences and all 4^k k-mers
bruteforceMedianString seqs k = map snd $ filter f $ zip totalDistances words
                          where totalDistances = map (total_dH seqs) words
                                words =  replicateM k dnaAlphabet
                                f (n, _) = n == minimum totalDistances
dnaAlphabet = "ACTG"

-- Hamming distance
dH seq1 seq2 = hamming seq1 seq2 0
    where hamming _ [] c = c
          hamming [] _ c = c 
          hamming (x:xs) (y:ys) c
              | x == y  = hamming xs ys c
              | otherwise = hamming xs ys $! (c + 1)

-- Total Hamming Distance
total_dH seqs k_mer = sum $ map min_dH seqs
                      where min_dH seq = minimum $ map (dH k_mer) $ k_mers1 (length k_mer) seq

-- All k-mers from sequence
-- For exapmle: 
-- k_mers 3 "AACCTT" 
-- ["AA","AC","CC","CT","TT"]
k_mers k seq = map (take k) $ take (n-k+1) $ tails seq
              where n = length seq

k_mers1 k seq = drop (k-1) $ tails' seq []
                where tails' []  acc = acc
                      tails' seq acc = tails' (tail seq) $ (take k seq):acc
                      n = length seq

-- Branch-and-Bound median string search.

data TreeItem = Item {k_mer :: String,
                      distance :: Int,
                      prefixDist :: Int}
                deriving (Show, Eq)

instance Ord TreeItem where
      compare i1 i2 = compare (distance i1) (distance i2)

-- Find motifs of length k in list of DNA sequences using branch-and-bound algorithm.
bbMedianString seqs k = improve k sTree
                        where sTree = searchTree seqs k

-- Improving search tree by pruning branches that begins in vertex 
-- with distance greater then or equal to current minimum distance.
improve k sTree | null $ availableVertices prunedTree k  = k_mer minItem
                | otherwise = improve k prunedTree
                  where minItem      = minimum $ map rootLabel nextVertices
                        prunedTree   = prune (>= minItem) sTree
                        nextVertices = availableVertices sTree k 

-- Generating serch tree. Node is represented as TreeItem data type.
searchTree seqs k = unfoldTree (nextLevel seqs k) $ Item "" 0 0

makeItem _ [] = []
makeItem prevItem (x:xs) = f x:makeItem (f x) xs 
                           where f x | dist < (distance prevItem) = Item x dist dist
                                     | otherwise = Item x (distance prevItem) dist
                                     where dist = total_dH sampleSeqs x

-- Consequentially adding all character from dnaAlphabets to prefix.
nextLevel seqs k vertex | level < k  = (vertex, map makeItem $ map (appendChar $ k_mer vertex) dnaAlphabet)
                        | otherwise  = (vertex, [])
                     where appendChar str c = str ++ [c]
                           makeItem x = Item x (dist $ splitAt middle x) (fixedDistance vertex)
                                  where dist (prefix, []) = total_dH seqs prefix 
                                        dist (prefix, suffix)
                                            | level == k-1 = total_dH seqs x                                            
                                            | otherwise = prefixDist vertex + (total_dH seqs suffix)
                                        fixedDistance vertex
                                            | level == middle = total_dH seqs $ k_mer vertex
                                            | otherwise       = prefixDist vertex                               
                           middle = k `div` 2
                           level = length $ k_mer vertex

-- Find vertices with potential solution.
availableVertices tree k | null nextToLastVertices = []
                         | otherwise = subForest $ head $ nextToLastVertices
                           where nextToLastVertices | null vertices  = [] 
                                                    | otherwise      =  dropWhile (null . subForest) $ head vertices 
                                                      where vertices = drop (k-1) (levels' tree)

sampleSeqs = ["TGACCGTGCCCTTGGA", "CCCTTGGAAGAAAAATGG", "AAACCTTGGACATGACT"]

main = do
    args <- getArgs
    case args of
        [fileName, k] -> do
           seqs <- readFasta fileName 
           putStrLn $ bbMedianString (map (toStr . seqdata) seqs) (read k)
        _  -> putStrLn "Wrong arguments count. Usage: MotifSearch <filename> <k-mer length>"    
    
    