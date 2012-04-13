import Data.Char (digitToInt)
import System.Random
import Control.Monad.State

asInt_fold [] = 0
asInt_fold ['-'] = 0
asInt_fold ('-':xs) = - asInt_fold xs
asInt_fold xs = foldl addDigit 0 xs
                where addDigit acc x = acc * 10 + digitToInt x

safeHead [] = Nothing
safeHead xs = Just (head xs)

myConcat [] = []
myConcat xs = foldr step [] xs
              where step xs ys = foldr (:) ys xs

myTakeWhile p (x:xs) = if p x then x:myTakeWhile p xs 
                       else []

myTakeWhile_foldr p xs = foldr step [] xs
                         where step x xs |p x = x:xs  
                                         |otherwise = []

data FastaFormat = Fasta {header, sequence :: String}
                   deriving (Show)

myReverse [] = []
myReverse (x:xs) = (myReverse xs) ++ [x]

myReverse1 xs = foldl (flip (:)) [] xs

type RandomState a = State StdGen a

getRandom :: Random a => RandomState a
getRandom =
  get >>= \gen ->
  let (val, gen') = random gen in
  put gen' >>
  return val

rnd :: StdGen -> (Int, StdGen)
rnd = random

getTwoRandoms :: Random a => RandomState (a, a)
getTwoRandoms = liftM2 (,) getRandom getRandom

getTwoIntRnd :: RandomState (Int,Int)
getTwoIntRnd = getTwoRandoms

randoms' :: StdGen -> [Int]
randoms' gen = let (value, newGen) = random gen in value:randoms' newGen