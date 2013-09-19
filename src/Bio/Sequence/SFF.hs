{- | Read and write the SFF file format used by
   Roche\/454 sequencing to store flowgram data.

   A flowgram is a series of values (intensities) representing homopolymer runs of
   A,G,C, and T in a fixed cycle, and usually displayed as a histogram.

   This file is based on information in the Roche FLX manual.  Among other sources for information about
   the format, are The Staden Package, which contains an io_lib with a C routine for parsing this format.
   According to comments in the sources, the io_lib implementation is based on a file
   called getsff.c, which I've been unable to track down.  Other software parsing SFFs 
   are QIIME, sff_extract, and Celera's sffToCa.

   It is believed that all values are stored big endian.
-}

module Bio.Sequence.SFF ( SFF(..), CommonHeader(..)
                        , ReadHeader(..), ReadBlock(..)
                        , readSFF, writeSFF, writeSFF', recoverSFF
                        -- , sffToSequence, rbToSequence
                        , trim, trimFromTo -- , trimKey
                        , baseToFlowPos, flowToBasePos
                        , trimFlows
                        , test, convert, flowgram
                        , masked_bases, cumulative_index
                        , packFlows, unpackFlows
                        , Flow, Qual, Index, SeqData, QualData
                        , ReadName (..), decodeReadName, encodeReadName
                        , putRB, getRB
                        ) where

import Bio.Core.Sequence
import Bio.Sequence.SFF_name

import Data.Int
import qualified Data.ByteString.Lazy as LB
import qualified Data.ByteString.Lazy.Char8 as LBC
import qualified Data.ByteString as B
import qualified Data.ByteString.Char8 as BC
import Data.ByteString (ByteString)
import Control.Monad (when,replicateM,replicateM_)

import Data.List (intersperse)
import Data.Binary
import Data.Binary.Get (getByteString,getLazyByteString,runGetState)
import qualified Data.Binary.Get as G
import Data.Binary.Put (putByteString,putLazyByteString)
import Data.Char (toUpper, toLower)
import Text.Printf (printf)
import System.IO

-- | The type of flowgram value
type Flow = Int16
type Index = Word8

-- Global variables holding static information
-- | An SFF file always start with this magic number.
magic :: Int32
magic = 0x2e736666

-- | Version is always 1.
versions :: [Int32]
versions = [1]

-- | Read an SFF file.
readSFF :: FilePath -> IO SFF
readSFF f = do
  file <- LB.readFile f
  let (header, remaining, _) = runGetState (get::Get CommonHeader) file 0
      blocks = getBlocks header (fromIntegral $ num_reads header) remaining
  return (SFF header blocks)

getBlocks :: CommonHeader -> Int -> LB.ByteString -> [ReadBlock]
getBlocks header n str
  | n == 0 = []
  | otherwise = case runGetState (getBlock $ fromIntegral $ flow_length $ header) str 0 of
    (block, remaining, _) -> block : getBlocks header (n-1) remaining

getBlock :: Int -> Get ReadBlock
getBlock flows = get >>= getRB flows

{-
-- | Extract the read without the initial (TCAG) key.
trimKey :: CommonHeader -> Sequence Nuc -> Maybe (Sequence Nuc)
trimKey ch (Seq n s q) = let (k,s2) = LB.splitAt (fromIntegral $ key_length ch) s
                          in if LBC.map toLower k==LBC.map toLower (LB.fromChunks [key ch]) 
                             then Just $ Seq n s2 (liftM (LB.drop (fromIntegral $ key_length ch)) q)
                             else Nothing -- error ("Couldn't match key in sequence "++LBC.unpack n++" ("++LBC.unpack k++" vs. "++BC.unpack (key ch)++")!")
-}

instance BioSeq ReadBlock where
  seqid rb = SeqLabel $ LB.fromChunks [read_name $ read_header rb]
  seqheader rb = hdr
    where n = read_name $ read_header rb
          hdr = case decodeReadName n of
            Just h -> SeqLabel $ LBC.unwords 
                      [LB.fromChunks [n], LBC.pack $ show $ date h
                      ,LBC.pack $ show $ time h, LBC.pack $ show $ region h
                      ,LBC.pack $ (show (x_loc h) ++ ":"++show (y_loc h))]
            Nothing -> seqid rb
  seqdata  rb = let h = read_header rb
                    (left,right) = (clip_qual_left h, clip_qual_right h)
                    (a,b) = LB.splitAt (fromIntegral right) $ unSD $ bases rb
                    (c,d) = LB.splitAt (fromIntegral left-1) a
                in SeqData $ LBC.concat [LBC.map toLower c, LBC.map toUpper d,LBC.map toLower b]
  seqlength rb = fromIntegral $ num_bases $ read_header rb

instance BioSeqQual ReadBlock where
  seqqual = quality

{- -- | Extract the sequences from an 'SFF' data structure.
sffToSequence :: SFF -> [Sequence Nuc]
sffToSequence (SFF _ rs) = map rbToSequence rs

-- | Extract the sequence information from a 'ReadBlock'.
rbToSequence :: ReadBlock -> Sequence Nuc
rbToSequence r = Seq (LB.fromChunks [read_name h ,BC.pack (" qclip: "++show left ++".."++show right)])
                     (seqdata r)
                     (Just $ seqqual r)
  where h = read_header r
        (left,right) = (clip_qual_left h, clip_qual_right h)
-}

-- | Trim a 'ReadBlock' limiting the number of flows.  If writing to
--   an SFF file, make sure you update the 'CommonHeader' accordingly.
--   See @examples/Flx.hs@ for how to use this.  
trimFlows :: Integral i => i -> ReadBlock -> ReadBlock
trimFlows l rb = rb { read_header = rh { num_bases = fromIntegral n
                                       , clip_qual_right = min cqr $ fromIntegral n
                                       }
                    , flow_data   = B.take (2*fromIntegral l) (flow_data rb)
                    , flow_index  = B.take n (flow_index rb)
                    , bases       = SeqData $ LB.take (fromIntegral n) (unSD $ bases rb)
                    , quality     = QualData $ LB.take (fromIntegral n) (unQD $ quality rb)
                    }
  where n = (flowToBasePos rb l)-1
        rh = read_header rb
        cqr = clip_qual_right rh

-- trimming the flowgram is necessary, but how to deal with the shift in flow
-- sequence - i.e. what to do when trimming "splits" a flow into trimmed/untrimmed bases?

-- | Trim a read to specific sequence position, inclusive bounds.
trimFromTo :: (Integral i) => i -> i -> ReadBlock -> ReadBlock
trimFromTo x r rd = let
  l = x-1
  trim_seq = LB.drop (fromIntegral l) . LB.take (fromIntegral r)
  trim_seq' = B.drop (fromIntegral l) . B.take (fromIntegral r)
  trim_flw = B.drop ((2*) $ fromIntegral $ baseToFlowPos rd l) . B.take ((2*) $ fromIntegral $ baseToFlowPos rd r)
  new_flw  = trim_flw (flow_data rd)
  padding = B.replicate (B.length (flow_data rd) - B.length new_flw) 0
  rh = read_header rd
  [r',l'] = map fromIntegral [r,l]
  rh' = rh { num_bases = fromIntegral (r'-l')
           , clip_qual_left = max 0 $ clip_qual_left rh-l'
           , clip_qual_right = min (clip_qual_right rh-l') (r'-l'+1)
           , clip_adapter_left = max 0 $ clip_adapter_left rh-l'
           , clip_adapter_right = min (clip_adapter_right rh-l') (r'-l'+1)
           }
  in rd { read_header = rh'
        , flow_data = B.concat [new_flw, padding]
        , flow_index = trim_seq' (flow_index rd)
        , bases = SeqData $ trim_seq $ unSD $ bases rd
        , quality = QualData $ trim_seq $ unQD $ quality rd
        }

-- | Trim a read according to clipping information
trim :: ReadBlock -> ReadBlock
trim rb = let rh = read_header rb in trimFromTo (clip_qual_left rh) (clip_qual_right rh) rb

-- | Convert a flow position to the corresponding sequence position
flowToBasePos :: Integral i => ReadBlock -> i -> Int
flowToBasePos rd fp = length $ takeWhile (<=fp) $ scanl (+) 0 $ map fromIntegral $ B.unpack $ flow_index rd

-- | Convert a sequence position to the corresponding flow position
baseToFlowPos :: Integral i => ReadBlock -> i -> Int
baseToFlowPos rd sp = sum $ map fromIntegral $ B.unpack $ B.take (fromIntegral sp) $ flow_index rd

-- | Read an SFF file, but be resilient against errors.
recoverSFF :: FilePath -> IO SFF
recoverSFF f = return . unRecovered . decode =<< LB.readFile f

-- | Write an 'SFF' to the specified file name
writeSFF :: FilePath -> SFF -> IO ()
writeSFF = encodeFile

-- | Write an 'SFF' to the specified file name, but go back and
--   update the read count.  Useful if you want to output a lazy
--   stream of 'ReadBlock's.  Returns the number of reads written.
writeSFF' :: FilePath -> SFF -> IO Int
writeSFF' f (SFF hs rs) = do
  h <- openFile f WriteMode
  LBC.hPut h $ encode hs
  c <- writeReads h (fromIntegral $ flow_length hs) rs
  hSeek h AbsoluteSeek 20
  LBC.hPut h $ encode c
  hClose h
  return $ fromIntegral c

-- | Write 'ReadBlock's to a file handle.
writeReads :: Handle -> Int -> [ReadBlock] -> IO Int32
writeReads _ _ [] = return 0
writeReads h i xs = go 0 xs
  where go c (r:rs) = do
          LBC.hPut h $ encode (RBI i r)
          let c' = c+1
          c' `seq` go c' rs
        go c [] = return c
        
data RBI = RBI Int ReadBlock

-- | Wrapper for ReadBlocks since they need additional information
instance Binary RBI where 
    put (RBI c r) = do
      putRB c r
    get = undefined
      
-- --------------------------------------------------
-- | test serialization by output'ing the header and first two reads 
--   in an SFF, and the same after a decode + encode cycle.
test :: FilePath -> IO ()
test file = do 
  (SFF h rs) <- readSFF file 
  let sff = (SFF h (take 2 rs))
  putStrLn $ show $ sff
  putStrLn ""
  putStrLn $ show $ (decode $ encode sff :: SFF)

-- --------------------------------------------------
-- | Convert a file by decoding it and re-encoding it
--   This will lose the index (which isn't really necessary)
convert :: FilePath -> IO ()
convert file = writeSFF (file++".out") =<< readSFF file

-- | Generalized function for padding
pad :: Integral a => a -> Put
pad x = replicateM_ (fromIntegral x) (put zero) where zero = 0 :: Word8 

-- | Generalized function to skip padding
skip :: Integral a => a -> Get ()
skip = G.skip . fromIntegral

-- | The data structure storing the contents of an SFF file (modulo the index)
data SFF = SFF !CommonHeader [ReadBlock]

instance Show SFF where 
    show (SFF h rs) = (show h ++ "Read Blocks:\n\n" ++ concatMap show rs)

instance Binary SFF where
    get = do
      -- Parse CommonHeader
      chead <- get
      -- Get the ReadBlocks
      rds <- replicateM (fromIntegral (num_reads chead))
                                   (do 
                                      rh <- get :: Get ReadHeader
                                      getRB (fromIntegral $ flow_length chead) rh
                                   )
      return (SFF chead rds)

    put (SFF hd rds) = do
      put hd
      mapM_ (put . RBI (fromIntegral $ flow_length hd)) rds

-- | Helper function for decoding a 'ReadBlock'.
{-# INLINE getRB #-}
getRB :: Int -> ReadHeader -> Get ReadBlock
getRB fl rh = do
  let nb = fromIntegral $ num_bases rh
      nb' = fromIntegral $ num_bases rh
  fg <- getByteString (2*fl)
  fi <- getByteString nb
  bs <- getLazyByteString nb'
  qty <- getLazyByteString nb'
  let l = (fl*2+nb*3) `mod` 8
  when (l > 0) (skip (8-l))
  return (ReadBlock rh fg fi (SeqData bs) (QualData qty))

-- | A ReadBlock can't be an instance of Binary directly, since it depends on
--   information from the CommonHeader.
putRB :: Int -> ReadBlock -> Put
putRB fl rb = do
  put (read_header rb)
  putByteString (flow_data rb)
  -- ensure that flowgram has correct lenght
  replicateM_ (2*fl-B.length (flow_data rb)) (put (0::Word8))
  putByteString (flow_index rb)
  putLazyByteString (unSD $ bases rb)
  putLazyByteString (unQD $ quality rb)
  let nb = fromIntegral $ num_bases $ read_header rb
      l = (fl*2+nb*3) `mod` 8
  when (l > 0) (pad (8-l))

-- | Unpack the flow_data field into a list of flow values
unpackFlows :: ByteString -> [Flow]
unpackFlows = dec . map fromIntegral . B.unpack 
    where dec (d1:d2:rest) = d1*256+d2 : dec rest
          dec [] = []
          dec _  = error "odd flowgram length?!"

-- | Pack a list of flows into the corresponding binary structure (the flow_data field)
packFlows :: [Flow] -> ByteString
packFlows = B.pack . map fromIntegral . merge 
  where merge (x:xs) = let (a,b) = x `divMod` 256 in a:b:merge xs
        merge [] = []

-- ----------------------------------------------------------
-- | SFF has a 31-byte common header
--
--   The format is open to having the index anywhere between reads,
--   we should really keep count and check for each read.  In practice, it
--   seems to be places after the reads.
--   
--   The following two fields are considered part of the header, but as
--   they are static, they are not part of the data structure
--
-- @        
--     magic   :: Word32   -- 0x2e736666, i.e. the string \".sff\"
--     version :: Word32   -- 0x00000001
-- @
data CommonHeader = CommonHeader {
          index_offset                            :: Int64    -- ^ Points to a text(?) section
        , index_length, num_reads                 :: Int32
        , key_length, flow_length                 :: Int16
        , flowgram_fmt                            :: Word8
        , flow, key                               :: ByteString 
        }

instance Show CommonHeader where
    show (CommonHeader io il nr kl fl fmt f k) =
        "Common Header:\n\n" ++ (unlines $ map ("    "++) 
                                 ["index_off:\t"++show io ++"\tindex_len:\t"++show il
                                 ,"num_reads:\t"++show nr
                                 ,"key_len:\t"  ++show kl ++ "\tflow_len:\t"++show fl
                                 ,"format\t:"   ++show fmt
                                 ,"flow\t:"     ++BC.unpack f
                                 ,"key\t:"      ++BC.unpack k
                                 , ""
                                 ])

instance Binary CommonHeader where
    get = do { m <- get ; when (m /= magic)   $ error (printf "Incorrect magic number - got %8x, expected %8x" m magic)
             ; v <- get ; when (not (v `elem` versions)) $ error (printf "Unexpected version - got %d, supported are: %s" v (unwords $ map show versions))
             ; io <- get ; ixl <- get ; nrd <- get
             ; chl <- get ; kl <- get ; fl <- get ; fmt <- get
             ; fw <- getByteString (fromIntegral fl)
             ; k  <- getByteString (fromIntegral kl)
             ; skip (chl-(31+fl+kl)) -- skip to boundary
             ; return (CommonHeader io ixl nrd kl fl fmt fw k)
             }

    put ch = let CommonHeader io il nr kl fl fmt f k = ch { index_offset = 0 } in
        do { let cl = 31+fl+kl
                 l = cl `mod` 8
                 padding = if l > 0 then 8-l else 0
           ; put magic; put (last versions); put io; put il; put nr; put (cl+padding); put kl; put fl; put fmt
           ; putByteString f; putByteString k
           ; pad padding -- skip to boundary
           }

-- ---------------------------------------------------------- 
-- | Each Read has a fixed read header, containing various information.
data ReadHeader = ReadHeader {
      name_length                           :: Int16
    , num_bases                             :: Int32
    , clip_qual_left, clip_qual_right
    , clip_adapter_left, clip_adapter_right :: Int16
    , read_name                             :: ByteString
}

instance Show ReadHeader where
    show (ReadHeader nl nb cql cqr cal car rn) =
        ("    Read Header:\n" ++) $ unlines $ map ("        "++) 
                    [ "name_len:\t"++show nl, "num_bases:\t"++show nb
                    , "clip_qual:\t"++show cql++"..."++show cqr
                    , "clip_adap:\t"++show cal++"..."++show car
                    , "read name:\t"++BC.unpack rn
                    , "" 
                    ]

instance Binary ReadHeader where
    get = do
      { rhl <- get; nl <- get; nb <- get
      ; cql <- get; cqr <- get ; cal <- get ; car <- get
      ; n <- getByteString (fromIntegral nl)
      ; skip (rhl - (16+ nl))
      ; return (ReadHeader nl nb cql cqr cal car n)
      }
    put (ReadHeader nl nb cql cqr cal car rn) = 
        do { let rl = 16+nl
                 l = rl `mod` 8
                 padding = if l > 0 then 8-l else 0
           ; put (rl+padding); put nl; put nb; put cql; put cqr; put cal; put car
           ; putByteString rn 
           ; pad padding
           }

-- ----------------------------------------------------------
-- | This contains the actual flowgram for a single read.
data ReadBlock = ReadBlock {
      read_header                :: ! ReadHeader
    -- The data block
    , flow_data                  :: ! ByteString -- nb! use unpackFlows for this
    , flow_index                 :: ! ByteString
    , bases                      :: ! SeqData
    , quality                    :: ! QualData
    }

-- | Helper function to access the flowgram
flowgram :: ReadBlock -> [Flow]
flowgram = unpackFlows . flow_data

-- | Extract the sequence with masked bases in lower case
masked_bases :: ReadBlock -> SeqData
masked_bases rb = let
  l = fromIntegral $ clip_qual_left $ read_header rb
  r = fromIntegral $ clip_qual_right $ read_header rb
  SeqData s = bases rb
  in SeqData $ LBC.concat [ LBC.map toLower $ LBC.take (l-1) s
                , LBC.take r (LBC.drop (l-1) s)
                , LBC.map toLower $ LBC.drop r s]

-- | Extract the index as absolute coordinates, not relative.
cumulative_index :: ReadBlock -> [Int]
cumulative_index = scanl1 (+) . map fromIntegral . B.unpack . flow_index

instance Show ReadBlock where
    show (ReadBlock h f i (SeqData b) (QualData q)) =
        show h ++ unlines (map ("     "++) 
            ["flowgram:\t"++show (unpackFlows f)
            , "index:\t"++(concat . intersperse " " . map show . B.unpack) i
            , "bases:\t"++LBC.unpack b
            , "quality:\t"++(concat . intersperse " " . map show . LB.unpack) q
            , ""
            ])

-- ------------------------------------------------------------
-- | RSFF wraps an SFF to provide an instance of Binary with some more error checking.
data RSFF = RSFF { unRecovered :: SFF }

instance Binary RSFF where 
    get = do
      -- Parse CommonHeader
      chead <- get
      -- Get the first read block
      r1 <- do rh <- get 
               getRB (fromIntegral $ flow_length chead) rh
      -- Get subsequent read blocks
      rds <- replicateM (fromIntegral (num_reads chead))
                                   (do rh <- getSaneHeader (take 4 $ BC.unpack $ read_name $ read_header r1)
                                       getRB (fromIntegral $ flow_length chead) rh)
      return (RSFF $ SFF chead (r1:rds))
    put = error "You should not serialize an RSFF"

-- | This allows us to decode the constant parts of the read header for verifying its correcness.
data PartialReadHeader = PartialReadHeader {
      _pread_header_lenght                    :: Int16 
    , _pname_length                           :: Int16
    , _pnum_bases                             :: Int32
    , _pclip_qual_left, _pclip_qual_right
    , _clip_adapter_left, _pclip_adapter_right :: Int16
    , _pread_name                              :: ByteString -- length four
}

instance Binary PartialReadHeader where
    get = do { rhl <- get; nl <- get; nb <- get; ql <- get; qr <- get; al <- get; ar <- get; rn <- getByteString 4 
             ; return (PartialReadHeader rhl nl nb ql qr al ar rn) }
    put = error "You should not serialize a PartialReadHeader"

-- | Ensure that the header we're decoding matches our expectations.
getSaneHeader :: String -> Get ReadHeader
getSaneHeader prefix = do
  buf <- getLazyByteString 20
  decodeSaneH prefix buf  

-- | Decode a 'ReadHeader', verifying that the data make sense.
decodeSaneH :: String -> LBC.ByteString -> Get ReadHeader
decodeSaneH prefix buf = do
  let PartialReadHeader rhl nl _nb _ql _qr _al _ar rn = decode buf
  if rhl >= 20 && nl > 0 && all id (zipWith (==) prefix (BC.unpack rn))
      then do buf2 <- getLazyByteString (fromIntegral rhl-20)
              return (decode $ LB.concat [buf,buf2])
      else do x <- getLazyByteString 1 -- error "skip one byte, try again"
              decodeSaneH prefix (LBC.concat [buf,x])

