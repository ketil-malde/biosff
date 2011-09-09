-- | This implements a number of filters used in the Titanium pipeline, 
--   based on published documentation.
module Bio.Sequence.SFF_filters where

import Bio.Sequence.SFF (ReadBlock(..), ReadHeader(..)
                        , flowToBasePos, flowgram, cumulative_index)

import Bio.Core.Sequence
import qualified Data.ByteString.Lazy as B
import qualified Data.ByteString as SB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List (tails)
import Data.Char (toUpper)

-- Ti uses a set of filters, described in the (something) manual.
-- (GS Run Processor Application, section 3.2.2++)

-- ** Discarding filters

-- | DiscardFilters determine whether a read is to be retained or discarded
type DiscardFilter = ReadBlock -> Bool -- True to retain, False to discard

-- | This filter discards empty sequences.
discard_empty :: DiscardFilter
discard_empty rb = num_bases (read_header rb) >= 5

-- | Discard sequences that don't have the given key tag (typically TCAG) at the start
--   of the read.
discard_key :: String -> DiscardFilter
discard_key key rb = (map toUpper key==) $ take (length key) $ BL.unpack $ unSD $ bases rb

-- | 3.2.2.1.2 The "dots" filter discards sequences where the last positive flow is 
--   before flow 84, and flows with >5% dots (i.e. three successive noise values) 
--   before the last postitive flow.  The percentage can be given as a parameter.
discard_dots :: Double -> DiscardFilter
discard_dots p rb = let dotcount = SB.length $ SB.filter (>3) $ flow_index rb
                    in fromIntegral dotcount / fromIntegral (BL.length $ unSD $ bases rb) < p
                       && last (cumulative_index rb) >= 84

-- | 3.2.2.1.3 The "mixed" filter discards sequences with more than 70% positive flows.  
--   Also, discard with <30% noise, >20% middle (0.45..0.75) or <30% positive.
discard_mixed :: DiscardFilter
discard_mixed rb = let fs = dropWhile (<50) . reverse . flowgram $ rb
                       fl = dlength fs
                   in and
                      [ (dlength (filter (>50) fs) / fl) < 0.7 -- 70% positive
                      , (dlength (filter (<45) fs) / fl) > 0.3 -- 30% noise
                      , (dlength (filter (>75) fs) / fl) > 0.3 -- 30% postivie
                      , (dlength (filter (\f -> f<=75 && f>=45) fs) / fl) < 0.2
                      ]

-- | Discard a read if the number of untrimmed flows is less than n (n=186 for Titanium)
discard_length :: Int -> DiscardFilter
discard_length n rb = length (flowgram rb) >= n

-- ** Trimming filters

-- | TrimFilters modify the read, typically trimming it for quality
type TrimFilter = ReadBlock -> ReadBlock

-- | 3.2.2.1.4 Signal intensity trim - trim back until <3% borderline flows (0.5..0.7).
--   Then trim borderline values or dots from the end (use a window).
trim_sigint :: TrimFilter
trim_sigint rb = clipSeq rb (sigint rb)

-- n counts the "bad" flow values, m counts flow position
sigint :: ReadBlock -> Int
sigint rb = let bs = drop 1 $ scanl (\(n,m,_) f -> if f >= 50 && f <= 70 then (n+1,m+1,f) else (n,m+1,f)) (0,0,0) $ flowgram rb 
                xs = dropWhile (\(_,_,f) -> f<=70) 
                     $ dropWhile (\(n,m,_)->(1000*n) `div` m > (30::Int)) 
                     $ reverse bs
            in case xs of []          -> error "no sequence left?"
                          ((_,m,_):_) -> flowToBasePos rb m

-- | 3.2.2.1.5 Primer filter 
-- This looks for the B-adaptor at the end of the read.  The 454 implementation isn't very
-- effective at finding mutated adaptors.
trim_primer :: String -> TrimFilter
trim_primer s rb = clipSeq rb (find_primer s rb)

find_primer :: String -> ReadBlock -> Int
find_primer s rb = go (num_bases (read_header rb) - 10)
  where go i | i <= 5    = fromIntegral (num_bases $ read_header rb)
             | match i   = fromIntegral i
             | otherwise = go (i-1)
        match j = s' `B.isPrefixOf` B.drop (fromIntegral j) (unSD $ bases rb)
        s' = BL.pack $ map toUpper $ take 14 s

-- 3.2.2.1.6 Trimback valley filter is ignored, we don't understand the description.

-- | 3.2.2.1.7 Quality score trimming trims using a 10-base window until a Q20 average is found.
trim_qual20 :: Int -> TrimFilter
trim_qual20 w rs = clipSeq rs $ qual20 w rs

qual20 :: Int -> ReadBlock -> Int
qual20 w rs = (fromIntegral $ num_bases $ read_header rs)
              - (length . takeWhile (<20) . map (avg . take w) . tails . reverse . B.unpack $ unQD $ quality rs)

-- ** Utility functions

-- | List length as a double (eliminates many instances of fromIntegral)
dlength :: [a] -> Double
dlength = fromIntegral . length

-- | Calculate average of a list
avg :: Integral a => [a] -> Double
avg xs = sum (map fromIntegral xs) / dlength xs

-- | Translate a number of flows to position in sequence, and update clipping data accordingly
clipFlows :: ReadBlock -> Int -> ReadBlock
clipFlows rb n = clipSeq rb (flowToBasePos rb n)

-- | Update clip_qual_right if more severe than previous value
clipSeq :: ReadBlock -> Int -> ReadBlock
clipSeq rb n' = let n = fromIntegral n' 
                    rh = read_header rb 
                in if clip_qual_right rh <= n then rb else rb { read_header = rh {clip_qual_right = n }}

-- ** Data

-- Celera docs, at http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=SffToCA

-- These are used for mate-pair libraries, should be located around the middle of the read:

flx_linker = "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC"  -- Celera
ti_linker  = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"  -- 20K cod jump

-- ti_linker and this: "AGCATATTGAAGCATATTACATACGATATGCTTCAATAATGC"
-- from "GS FLX Titanium 3 kb Span Paired End Library Preparation Method Manual April 2009"
-- ftp://ftp.genome.ou.edu/pub/for_broe/titanium/

-- These are used at the end of RNA (cDNA) sequences, after the poly-A tail:

rna_adapter   = "ggcgggcgatgtctcgtctgagcgggctggcaaggc" -- cod transcripts?
rna_adapter2  = "ttcgcagtgagtgacaggctagtagctgagcgggctggcaaggc"  -- Cod_c.sff
rna_adapter3  = "gacggggcggatgtctcgtctgagcgggcgtggcaaggc"       -- COD1.sff

-- These are used at the end of DNA sequencing reads:

rapid_adapter = "agtcgtggaggcaaggcacacagggatagg"  -- sea louse reads, key GACT
ti_adapter_b  = "ctgagactgccaaggcacacagggggatagg"  -- sea bass and l.s.Ca
                 


