{-# LANGUAGE DeriveDataTypeable #-}

module Options where

import System.Console.CmdArgs
import Control.Monad (when)
import Data.Maybe (isJust, isNothing)

data Opts = Opts 
            { trimKey :: Bool
            , trimAdapter :: Bool
            , trim    :: Bool
            , filterEmpty :: Bool
            , plot    :: Maybe String
            , summarize :: Maybe FilePath
            , filters :: Maybe FilePath
            , info    :: Maybe FilePath
            , fasta   :: Maybe FilePath
            , fqual  :: Maybe FilePath
            , fastq   :: Maybe FilePath
            , flowgram :: Maybe FilePath
            , histogram :: Maybe FilePath
            , histpos :: Maybe FilePath
            , inputs :: [FilePath]
            , text   :: Maybe FilePath
            } deriving (Data,Typeable, Show, Eq)

optdef :: Ann
optdef = opt ("-"::String)

opts :: Opts
opts = Opts
  { trimKey = False &= help "Trim only the TCAG key sequence"
  , trim    = False &= help "Trim quality using clipping information"     &= name "t"
  , trimAdapter = False &= help "Trim quality using adapter information"     &= name "a"
  , filterEmpty = False &= help "Remove reads that are empty after trimming" &= name "E"
  , summarize = def   &= help "Output per sequence summary information"   &= typFile &= optdef
  , filters   = def   &= help "Output filtering information"              &= typFile &= optdef
  , info    = def   &= help "Output brief overview of the contents"       &= typFile &= optdef
  , fasta   = def   &= help "Output FASTA-formatted sequences"            &= typFile &= name "f" &= optdef
  , fqual    = def   &= help "Output phred qualities"                      &= typFile &= name "q" &= optdef
  , fastq = def   &= help "Output FastQ-formatted sequence and Sanger quality" &= typFile &= name "Q" &= optdef
  , flowgram = def  &= help "Output flowgram information in tabular form" &= typFile &= name "F" &= optdef
  , histogram = def &= help "Output histogram of flow values by nucleotide" &= typFile &= name "h" &= optdef
  , histpos   = def &= help "Output histogram of flow values by flow cycle" &= typFile &= name "H" &= optdef
  , plot      = Nothing &= help "Output gnuplot script for visualization" &= opt ""
  , text      = def &= help "Output SFF information as text (default)"    &= typFile &= name "T" &= optdef
  , inputs  = def &= args &= typFile
  } 
  &= summary "flower (biosff v0.3.7.1) - Extract information from SFF files" 
  &= program "flower"

getArgs :: IO Opts
getArgs = do
  o <- cmdArgs opts 
  -- print o
  let outs = filter isJust $ map ($o) [summarize,filters,info,fasta,fqual,fastq,flowgram,histogram,histpos,text]
  when ((length $ filter (==Just "-") $ outs) > 1) $ error "If you specify more than one output format, you need to specify output files"
  when (isJust (plot o) && isNothing (histogram o)) $ error "Gnuplot output only supported for histograms"
  let o' = if null outs then o { text = Just "-" } else o
  return o'
