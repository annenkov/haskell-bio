-- Tree utils as in "Modular Lazy Search for Constraint Satisfaction Problems" by  Thomas Nordin and Andrew Tolmach
module TreeUtils where
import Data.Tree

type Transform a b = Tree a -> Tree b

-- Like levels but returns [[Tree a]]
--levels' :: Tree a -> [[a]]
levels' t =  takeWhile (not . null) $
	     iterate (concatMap subForest) [t]

mapTree :: (a -> b) -> Transform a b
mapTree f (Node a cs) = Node (f a) (map (mapTree f) cs)

foldTree :: (a -> [b] -> b) -> Tree a -> b
foldTree f (Node a cs) = f a (map (foldTree f) cs)

filterTree :: (a -> Bool) -> Transform a a
filterTree p = foldTree f
               where f a cs = Node a (filter (p . rootLabel) cs)

prune :: (a -> Bool) -> Transform a a
prune p = filterTree (not . p)