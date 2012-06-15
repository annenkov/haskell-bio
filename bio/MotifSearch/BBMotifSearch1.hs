module BBMotifSearch1 where

import Control.Monad
import Data.List
import TreeUtils
import Data.Tree
import System (getArgs)
import Bio.Sequence.Fasta
import Bio.Sequence.SeqData
import BruteforceMotifSearch (total_dH, dnaAlphabet)

-- Improved Branch-and-Bound median string search.

data TreeItem = Item {k_mer :: String,
                      distance :: Int,
                      distances :: [Int]}
                deriving (Show, Eq)

instance Ord TreeItem where
      compare i1 i2 = compare (distance i1) (distance i2)

-- Find motifs of length k in list of DNA sequences using improved branch-and-bound algorithm.
bbMedianString seqs k = improve k sTree
                        where sTree = searchTree seqs k

-- Improving search tree by pruning branches that begins in vertex 
-- with distance greater then or equal to current minimum distance.
improve k sTree | null $ availableVertices prunedTree k  = k_mer minItem
                | otherwise = improve k prunedTree
                  where minItem      = minimum $ map rootLabel nextVertices
                        prunedTree   = prune (>= minItem) sTree
                        nextVertices = availableVertices sTree k 

-- Generating search tree. Node is represented as TreeItem data type.
searchTree seqs k = unfoldTree (nextLevel seqs k) $ Item "" 0 []

upperLevel seqs k_mers@(k_mer:ks) Item {distances=ds} = upperLevel' (Item k_mer dist (dist:ds)) k_mers
                         where upperLevel' _ [] = []
                               upperLevel' Item {distances=dsts@(d:ds1)} (x:xs) = f x:upperLevel' (f x) xs
                                   where f x | dist < d   = Item x dist (dist:ds1)
                                             | otherwise  = Item x dist dsts
                                         dist = total_dH seqs x
                               dist = total_dH seqs k_mer


-- Consequentially adding all character from dnaAlphabets to prefix.
nextLevel seqs k vertex | level == k `div` 2 - 1 && even k = (vertex, map makeItem' $ map (appendChar $ k_mer vertex) dnaAlphabet)
                        | level < k `div` 2 = (vertex, upperLevel seqs (map (appendChar $ k_mer vertex) dnaAlphabet) vertex)
                        | level < k  = (vertex, map makeItem $ map (appendChar $ k_mer vertex) dnaAlphabet)
                        | otherwise  = (vertex, [])
                          where appendChar str c = str ++ [c]
--                                makeUpperLevel (x:xs) = firstVertex:upperLevel seqs firstVertex xs
--                                                        where firstVertex  = Item x dist (dist:(distances vertex))
--                                                              dist = (total_dH seqs x)
                                makeItem x | level < k-1 = Item x (dist x) (tail $ distances vertex)
                                           | otherwise = Item x (total_dH seqs x) []
                                             where dist x = (total_dH seqs x) + (head $ distances vertex)
                                makeItem' x = Item x ((total_dH seqs x) * 2) (distances vertex)
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
    
    