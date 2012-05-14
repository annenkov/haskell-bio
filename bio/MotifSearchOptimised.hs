module MotifSearch where

import Control.Monad
import Data.List
import TreeUtils
import Data.Tree
import qualified Data.ByteString.Char8 as BString
import System (getArgs)
import Bio.Sequence.Fasta
import Bio.Sequence.SeqData

dnaAlphabet = "ACTG"

-- Hamming distance
dH seq1 seq2 = sum $ BString.zipWith dist seq1 seq2
                where dist a b | a==b      = 0
                               | otherwise = 1

-- Total Hamming Distance
total_dH seqs k_mer = sum $ map min_dH seqs
                      where min_dH seq = minimum $ map (dH k_mer) $ k_mers (BString.length k_mer) seq

-- All k-mers from sequence
-- For exapmle: 
-- k_mers 3 "AACCTT" 
-- ["AA","AC","CC","CT","TT"]
k_mers k seq = map (BString.take k) $ take (n-k+1) $ BString.tails seq
               where n = BString.length seq

-- Branch-and-Bound median string search.

data TreeItem = Item {k_mer :: BString.ByteString,
                      distance :: Int}
                deriving (Show, Eq)

instance Ord TreeItem where
      compare i1 i2 = compare (distance i1) (distance i2)

-- Find motifs of length k in list of DNA sequences using branch-and-bound algorithm.
bbMedianString :: [BString.ByteString] -> Int -> BString.ByteString
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
searchTree seqs k = mapTree f $ unfoldTree (nextLevel k) BString.empty
                    where f x = Item x (total_dH seqs x)

-- Consequentially adding all character from dnaAlphabets to prefix.
nextLevel k prefix | length prefix < k      = (prefix, map (appendChar prefix) dnaAlphabet)
                   | otherwise              = (prefix, [])
                     where appendChar str c = str `BString.append` (BString.singleton c)

-- Find vertices with potential solution.
availableVertices tree k | null nextToLastVertices = []
                         | otherwise = subForest $ head $ nextToLastVertices
                           where nextToLastVertices | null vertices  = [] 
                                                    | otherwise      =  dropWhile (null . subForest) $ head vertices 
                                                      where vertices = drop (k-1) (levels' tree)

sampleSeqs = ["TGACCGTGCCCTTGGA", "CCTTGGAAGAAAAAATGG", "AAACCTTGGACATGACT"]

--main = do
--    args <- getArgs
--    case args of
 --       [fileName, k] -> do
 --          seqs <- readFasta fileName 
  --         putStrLn $ show $ bbMedianString (map seqdata seqs) (read k)
  --      _  -> putStrLn "Wrong arguments count. Usage: MotifSearch <filename> <k-mer length>"    
