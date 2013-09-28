{-| FRecover reads a possibly corrupt SFF file, and
  tries to extract as much information as possible from it,
  generating a new SFF file in the the process
-}
module Main where
import Bio.Sequence.SFF
import System.Environment (getArgs)

main :: IO ()
main = do 
  xs <- getArgs
  if null xs 
    then error "frecover (biosff 0.3.7.1) - try to repair broken SFF files.\nUsage: frecover file.sff [file.sff]" 
    else do mapM_ recoverFile xs
  where recoverFile f = writeSFF (f++"_recovered") =<< recoverSFF f
  -- perhaps we should use writeSFF' instead, since it goes back
  -- to update the header with number of reads written?
