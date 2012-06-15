-- Tree utils as in "Modular Lazy Search for Constraint Satisfaction Problems" by  Thomas Nordin and Andrew Tolmach
module TreeUtils where
import Data.Tree

type Transform a b = Tree a -> Tree b

levels' :: Tree a -> [Forest a]
levels' t =  takeWhile (not . null) $
	     iterate (concatMap subForest) [t]

foldTree :: (a -> [b] -> b) -> Tree a -> b
foldTree f (Node a cs) = f a (map (foldTree f) cs)

filterTree :: (a -> Bool) -> Transform a a
filterTree p = foldTree f
               where f a cs = Node a (filter (p . rootLabel) cs)

prune :: (a -> Bool) -> Transform a a
prune p = filterTree (not . p)