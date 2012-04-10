import Data.Char (digitToInt)

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