
module BBMotifSearch where

import Control.Monad
import Control.Monad.State
import Data.List
import TreeUtils
import Data.Tree
import System.Environment (getArgs)
import Bio.Sequence.Fasta
import Bio.Sequence.SeqData
import BruteforceMotifSearch (total_dH, dnaAlphabet, bruteforceMedianString)

-- Branch-and-Bound median string search.

data TreeItem = Item {
          k_mer    :: String
        , distance :: Int
        } deriving (Show, Eq)

instance Ord TreeItem where
      compare i1 i2 = compare (distance i1) (distance i2)

-- Find motifs of length k in list of DNA sequences using branch-and-bound algorithm.
bbMedianString seqs k = improve k sTree
                        where sTree = searchTree seqs k

-- Improving search tree by pruning branches that begins in vertex 
-- with distance greater then or equal to current minimum distance.
improve k sTree | null $ availableVertices prunedTree k = k_mer minItem
                | otherwise = improve k prunedTree
                  where minItem      = minimum $ map rootLabel nextVertices
                        prunedTree   = prune (>= minItem) sTree
                        nextVertices = availableVertices sTree k 
-- last $ filter ((== k) . length) $ map k_mer 
preorderSearch seqs k = last $ filter ((== k) . length) $ map k_mer $ evalState (preorder k sTree) (maxBound :: Int)
                        where sTree = searchTree seqs k
                              lastLevel x = (length $ k_mer x) == k

preorder k (Node val subForest) = do dist <- get
                                     lst <- preorder' k subForest                                     
                                     return $ [val] ++ lst

preorder' :: Int -> [Tree TreeItem] -> State Int [TreeItem]
preorder' k forest = do                                               
                       dist <- get
                       tmp <- filterM distanceLT forest
                       lst <- mapM (newDistance k) tmp
                       return $ foldl' (++) [] lst
            where distanceLT x = do { dist <- get; return $ distance (rootLabel x) < dist }
--                  getDist x = do {dist <- get; return (fst x, dist)}

isLastLevel _ [] = False
isLastLevel k (t:ts) = (length $ k_mer $ rootLabel t) == k
               
newDistance :: Int -> Tree TreeItem -> State Int [TreeItem]
newDistance k node@(Node item forest) = do
                                           dist <- get
                                           let minItem = minimum $ map rootLabel forest
                                           if (length $ k_mer item) == k-1
                                               then 
                                                  if (distance minItem) < dist 
                                                  then do {put $ distance minItem; return [rootLabel node, minItem]}
                                                  else return [minItem]
                                               else preorder k node

-- Generating search tree; rootLabel of Node is represented as TreeItem data type.
searchTree seqs k = fmap f $ unfoldTree (nextLevel k) ""
                    where f x = Item x (total_dH seqs x)

-- Consequentially adding all character from dnaAlphabet to prefix.
nextLevel k prefix | length prefix < k      = (prefix, map (appendChar prefix) dnaAlphabet)
                   | otherwise              = (prefix, [])
                     where appendChar str c = str ++ [c]

-- Find vertices with potential solution.
availableVertices tree k
    | null lastLevels || null nextToLastVertices = []
    | otherwise = subForest $ head nextToLastVertices
    where
      lastLevels = drop (k-1) (levels' tree)
      nextToLastVertices = dropWhile (null . subForest) $ head lastLevels

sampleSeqs = ["TGACCGTGCCCTTGGA", "CCCTTGGAAGAAAAATGG", "AAACCTTGGACATGACT"]

main = do
    args <- getArgs
    case args of
        [fileName, k] -> do
           seqs <- readFasta fileName 
--           putStrLn $  bbMedianString (map (toStr . seqdata) seqs) (read k)
           putStrLn $ show $ bruteforceMedianString (map (toStr . seqdata) seqs) (read k)
        _  -> putStrLn "Wrong arguments count. Usage: MotifSearch <filename> <k-mer length>"
