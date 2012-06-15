module MotifTests where 

import Control.Monad.Random
import Test.HUnit
import BruteforceMotifSearch (bruteforceMedianString, dnaAlphabet)
import BBMotifSearch (bbMedianString)
import System (getArgs)

-- Random sequences for testing --
rndNucleotide :: Rand StdGen Char
rndNucleotide = do i <- getRandomR (0,3)
                   return $ dnaAlphabet!!i

rndDNA n = sequence $ replicate n rndNucleotide

-- generate m DNA sequences of length n
rndSeqs m n = sequence $ replicate m $ rndDNA n

implantMotif k_mer seq = do dna <- seq
                            i <- getRandomR (0,length dna - length k_mer - 1)
                            return $ take i dna ++ k_mer ++ drop i dna

makeTestSeqs m n k_mer =  sequence $ replicate m $ implantMotif k_mer $ rndDNA n


-- Tests --

bruteforceMedianStringTest = TestCase (assertEqual "Finding 4-mer:" [k_mer] (bruteforceMedianString testData 4))
        where k_mer = "AAAA"
              testData = evalRand (makeTestSeqs 5 20 k_mer) (mkStdGen 1)

brunchAndBoundMedianStringTest = TestCase (assertEqual "Finding 4-mer:" k_mer (bbMedianString testData 4))
        where k_mer = "AAAA"
              testData = evalRand (makeTestSeqs 5 20 k_mer) (mkStdGen 1)


-- run tests with "runTestTT tests"
tests = TestList [bruteforceMedianStringTest, brunchAndBoundMedianStringTest]

main = do
    [n, m, k_mer] <- getArgs
    let testData = evalRand (makeTestSeqs (read m) (read n) k_mer) (mkStdGen 2)
    putStrLn $ bbMedianString testData (length k_mer)