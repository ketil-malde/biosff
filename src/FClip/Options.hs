{-# LANGUAGE DeriveDataTypeable #-}

module Options where

import System.Console.CmdArgs

data Opts = Opts 
            { trimKey :: Bool
            , trimAdapter :: Bool
            , trim    :: Bool
            , filterEmpty :: Bool
            , output :: FilePath
            , input  :: FilePath
            } deriving (Data,Typeable, Show, Eq)

opts :: Opts
opts = Opts
  { trimKey = False &= help "Trim only the TCAG key sequence"             &= name "k"
  , trim    = False &= help "Trim quality using clipping information"     &= name "t"
  , trimAdapter = False &= help "Trim quality using adapter information"  &= name "a"
  , filterEmpty = False &= help "Filter out reads that are empty after trimming" &= name "E"
  , output = def &= help "SFF file to write" &= name "o"
  , input  = def &= args &= typFile
  } 
  &= summary "fclip (biosff v0.3.7) - Extract information from SFF files" 
  &= program "fclip"

getArgs :: IO Opts
getArgs = do
  o <- cmdArgs opts 
  -- todo: option checking
  return o

            
            