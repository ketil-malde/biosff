{-| FRecover reads a possibly corrupt SFF file, and
  tries to extract as much information as possible from it,
  generating a new SFF file in the the process
-}
module Main where
import Bio.Sequence.SFF
import System.Environment (getArgs)

main :: IO ()
main = mapM_ recoverFile =<< getArgs
  where recoverFile f = writeSFF (f++"_recovered") =<< recoverSFF f
  -- perhaps we should use writeSFF' instead, since it goes back
  -- to update the header with number of reads written?
