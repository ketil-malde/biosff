{-# LANGUAGE DeriveDataTypeable #-}

module Options where

import System.Console.CmdArgs
import Control.Monad (when)
import Data.Maybe (isJust)

data Opts = Opts 
            { trimKey :: Bool
            , trim    :: Bool
            , summarize :: Maybe FilePath
            , filters :: Maybe FilePath
            , info    :: Maybe FilePath
            , fasta   :: Maybe FilePath
            , fqual  :: Maybe FilePath
            , fastq   :: Maybe FilePath
            , illumina :: Maybe FilePath
            , flowgram :: Maybe FilePath
            , histogram :: Maybe FilePath
            , inputs :: [FilePath]
            , text   :: Maybe FilePath
            } deriving (Data,Typeable, Show, Eq)

optdef :: Ann
optdef = opt ("-"::String)

opts :: Opts
opts = Opts
  { trimKey = False &= help "Trim only the TCAG key sequence"
  , trim    = False &= help "Trim quality using clipping information"     &= name "t"
  , summarize = def   &= help "Output per sequence summary information"   &= typFile &= optdef
  , filters   = def   &= help "Output filtering information"              &= typFile &= optdef
  , info    = def   &= help "Output brief overview of the contents"       &= typFile &= optdef
  , fasta   = def   &= help "Output FASTA-formatted sequences"            &= typFile &= name "f" &= optdef
  , fqual    = def   &= help "Output phred qualities"                      &= typFile &= name "q" &= optdef
  , fastq = def   &= help "Output FastQ-formatted sequence and Sanger quality" &= typFile &= name "Q" &= optdef
  , illumina = def   &= help "Output FastQ-formatted sequence and Illumina quality" &= typFile &= name "I" &= optdef
  , flowgram = def  &= help "Output flowgram information in tabular form" &= typFile &= name "F" &= optdef
  , histogram = def &= help "Output histogram of flow values"             &= typFile &= name "h" &= optdef
  , text      = def &= help "Output SFF information as text (default)"    &= typFile &= name "T" &= optdef
  , inputs  = def &= args &= typFile
  } 
  &= summary "flower v0.7 - Extract information from SFF files" 
  &= program "flower"

getArgs :: IO Opts
getArgs = do
  o <- cmdArgs opts 
  print o
  let outs = filter isJust $ map ($o) [summarize,filters,info,fasta,fqual,fastq,illumina,flowgram,histogram,text]
  when ((length $ filter (==Just "-") $ outs) > 1) $ error "If you specify more than one output format, you need to specify output files"
  let o' = if null outs then o { text = Just "-" } else o
  return o'
