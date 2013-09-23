module Main where

import Bio.Sequence.SFF
import qualified Options as O
import Control.Monad (when)

main :: IO ()
main = do
  o <- O.getArgs
  when (null $ O.input o) $ error "Please provide an input file - or use --help for more information."
  when (null $ O.output o) $ error "Please provide an output file - or use --help for more information."  
  SFF h rs <- readSFF $ O.input o
  let t = buildTrim o
      f = if O.filterEmpty o then filter (\r -> (num_bases (read_header r)) > 0) else id
  _ <- writeSFF' (O.output o) (SFF h $ f $ map t rs) 
  return ()

buildTrim :: O.Opts -> (ReadBlock -> ReadBlock)
buildTrim o = if O.trimKey o then trimKey else id
              . if O.trim o then trim else id
              . if O.trimAdapter o then trimAdapter else id
