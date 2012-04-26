import Data.Tree

twoBranches x = (show x, [x,x])

dnaAlphabet = "ACTG"

nextLevel k prefix | length prefix < k      = (prefix, map (appendChar prefix) dnaAlphabet)
                   | otherwise              = (prefix, [])
                     where appendChar str c = str ++ [c]

searchTree k = unfoldTree (nextLevel k) ""
