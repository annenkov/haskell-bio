import Data.MultiSet (MultiSet, findMax, isSubsetOf, difference, null)
import qualified Data.MultiSet as MultiSet
import Data.List

makeDistances x ys = map (\y -> abs (x-y)) ys

allDistances xs = concatMap f (tails xs)
                  where f [] = []
                        f (x:xs) = makeDistances x xs

allDistances1 [] acc     = acc
allDistances1 (x:xs) acc = allDistances1 xs (acc ++ makeDistances x xs)
                           where makeDistances x ys = map (\y -> abs (x-y)) ys

--partialDigest :: MultiSet Int -> [Int]
partialDigest distances =  place (delete width distances) [0,width] width
                           where width = (maximum distances)

place [] acc _ = Just $ sort acc
place distances acc width  = if checkFit y1 acc distances
                             then place (newDistances y1 distances acc) (y1:acc) width
                             else if checkFit y2 acc distances
                                  then place (newDistances y2 distances acc) (y2:acc) width
                                  else Nothing
                             where newDistances x distances acc = MultiSet.toList (difference (MultiSet.fromList distances)
                                                                                              (MultiSet.fromList (makeDistances x acc)))
                                   (y1,y2) = alternatives distances width
-- MultiSet.fromList [1,2,3]

newDistances x distances acc = MultiSet.toList (difference (MultiSet.fromList distances)
                                                           (MultiSet.fromList (makeDistances x acc)))

alternatives xs width = (maximum xs, width - (maximum xs))

checkFit newPoint points distances = MultiSet.fromList (makeDistances newPoint points) `isSubsetOf`  MultiSet.fromList distances

-- nextStep distances points width = 