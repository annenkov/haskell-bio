import Data.Tree

dnaAlphabet = "ACTG"

nextLevel k prefix | length prefix < k      = (prefix, map (appendChar prefix) dnaAlphabet)
                   | otherwise              = (prefix, [])
                     where appendChar str c = str ++ [c]

searchTree k = unfoldTree (nextLevel k) ""

type Transform a b = Tree a -> Tree b
label :: Tree a -> a
label (Node lab _) = lab

foldTree :: (a -> [b] -> b) -> Tree a -> b
foldTree f (Node a cs) = f a (map (foldTree f) cs)

filterTree :: (a -> Bool) -> Transform a a
filterTree p = foldTree f
               where f a cs = Node a (filter (p . label) cs)

prune :: (a -> Bool) -> Transform a a
prune p = filterTree (not . p)