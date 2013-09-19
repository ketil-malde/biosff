-- FlowEr - FLOWgram ExtractoR
module Main (main) where

import Bio.Sequence.SFF
import Bio.Sequence.SFF_filters
import Bio.Core hiding (toText)

import Print as P
import Text.Printf

import System.IO (stdout, Handle, openFile, IOMode(..), hClose, hPutStrLn)

import Numeric (showFFloat)
import Data.Char (toLower)
import Data.List (intersperse)
import Data.ByteString.Char8 (unpack,ByteString)
import qualified Data.ByteString.Char8 as B
import qualified Data.ByteString as B1
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.ByteString.Lazy as L1

import Data.Array.Unboxed
import Data.Array.Unsafe -- unsafeFreeze is deprecated in D.A.ST
import Data.Array.ST hiding (unsafeFreeze)
import Control.Monad.ST
import Control.Monad.State

import Metrics
import qualified Options as O
import Options (Opts)
import Fork

main :: IO ()
main = do
  opts <- O.getArgs
  when (null $ O.inputs opts) $ error "Please provide an input file - or use --help for more information."
  forkAndWait $ buildActions opts

type Action  = IO ()
type Trimmer = ReadBlock -> ReadBlock

buildActions :: Opts -> [Action]
buildActions o = let
    inp = mapM readSFF (O.inputs o)
    tr = map (mkTrimmer o)
    ch (SFF h _) = h
    rs (SFF _ r) = r
    in snd $ flip runState [] $ do
      on (O.info o)      (\h -> mapM_ (hPutStrLn h . getHeader . ch) =<< inp)
      on (O.fasta o)     (\h -> mapM_ (L1.hPut h . L1.concat . map toFasta . tr . rs) =<< inp)
      on (O.fqual o)     (\h -> mapM_ (L1.hPut h . L1.concat . map toFastaQual . tr . rs) =<< inp)
      on (O.text o)      (\h -> mapM_ (hPutStrLn h . dumpText . tr . rs) =<< inp)
      on (O.fastq o)     (\h -> mapM_ (L1.hPut h . L1.concat . map toFastQ . tr . rs) =<< inp)
      on (O.summarize o) (\h -> mapM_ (L1.hPut h . summarize . tr . rs) =<< inp)  -- should we trim?
      on (O.filters o)   (\h -> mapM_ (L1.hPut h . sum_filters . rs) =<< inp)
      on (O.histogram o) (\h -> mapM_ (\(SFF c r) -> hPutStrLn h . (case O.plot o of Just str -> showPlot str 1000; Nothing -> showHist 9999)  ["A","C","G","T"] . histogram (B.unpack $ flow c) . map flowgram . tr $ r) =<< inp)
      on (O.histpos o)   (\h -> mapM_ (\(SFF _ r) -> hPutStrLn h . showHist 549 [] . histpos 549 . tr $ r) =<< inp)
      on (O.flowgram o)  (\h -> mapM_ (\(SFF c r) -> L1.hPut h . L1.fromChunks . intersperse (B.pack "\n") . concatMap (showread c) $ r) =<< inp)

on :: Maybe FilePath -> (Handle -> Action) -> State [Action] ()
on Nothing _    = return ()
on (Just f) act = modify $ (:) $ case f of 
  "-" -> act stdout
  _   -> do h <- openFile f WriteMode
            act h
            hClose h

mkTrimmer :: Opts -> Trimmer
mkTrimmer o = case (O.trimKey o, O.trim o, O.trimAdapter o) of
        (True,False,False) -> \r -> trimFromTo 5 (num_bases $ read_header r) r
        (False,True,False) -> trim
        (False,False,True) -> trimAdapter        
        (False,False,False) -> id
        _ -> error "Please specify only one of --trim, --trimAdapter, and --trimkey"

trimAdapter :: Trimmer
trimAdapter r = trimFromTo (clip_adapter_left rh) (if car == 0 then fromIntegral (num_bases rh) else car) r
  where rh = read_header r
        car = clip_adapter_right rh
-- ------------------------------------------------------------
-- No option - dump as text format
-- ------------------------------------------------------------
dumpText :: [ReadBlock] -> String
dumpText rs = concat . map toText $ rs
  where toText :: ReadBlock -> String
        toText r = concat [ gt, B.unpack (read_name rh), nl
                          , maybe "" ((\s->info++s++nl) . formatRN) $ decodeReadName (read_name rh)
                          , let (lf,rt) = (clip_adapter_left rh, clip_adapter_right rh) 
                            in if lf /= 0 || rt /= 0 then adapter ++ show lf ++ sp++ show rt ++ nl else ""
                          , clip,     show (clip_qual_left rh), sp, show (clip_qual_right rh), nl
                          , flows,    B.unpack $ B.unwords $ map fi $ flowgram r, nl
                          , idx,      unwords $ map show $ cumulative_index' r, nl
                          , base,     L.unpack (unSD $ seqdata r), nl
                          , qual,     unwords $ map show $ L1.unpack (unQD $ quality r), nl
                             ]
          where rh = read_header r
                gt = ">"
                nl = "\n"
                sp = " "
                info     = "  Info: \t"
                clip     = "  Clip: \t"
                adapter  = "  Adap: \t"
                flows    = "  Flows:\t"
                idx      = "  Index:\t"
                base     = "  Bases:\t"
                qual     = "  Quals:\t"
                formatRN (ReadName (yr,mo,dy) (h,m,s) r' x y) = 
                  printf "%4d-%02d-%02d %02d:%02d:%02d R%d (%d,%d)" yr mo dy h m s r' x y

cumulative_index' :: ReadBlock -> [Int]
cumulative_index' = scanl1 (+) . map fromIntegral . B1.unpack . flow_index

-- ------------------------------------------------------------
-- The -i option: Print header info
-- ------------------------------------------------------------
getHeader :: CommonHeader -> String
getHeader h = unlines ["Index:    \t" ++ show (index_offset h,index_length h)
                      ,"Num_reads:\t" ++ show (num_reads h)
                      ,"Num_flows:\t" ++ show (flow_length h)
                      ,"Key:      \t" ++ unpack (key h)
                      ]

-- ----------------------------------------------------------
-- The -s option: Summarize each read on one line
-- ----------------------------------------------------------

-- | Summarize each read on one line of output
summarize :: [ReadBlock] -> L.ByteString
summarize rs = do
  L.concat [ L.pack "# name........\tdate......\ttime....\treg\ttrim_l\ttrim_r\tx_loc\ty_loc\tlen\tK2\ttrimK2\tncount\tavgQ\ttravgQ\n"
           , P.toLazyByteString . mconcat . map sum1 $ rs]

-- todo: date and time are usually constants!
sum1 :: ReadBlock -> Builder
sum1 r = let rh = read_header r
             nb = num_bases rh
             h = read_name rh
             tr = trim r
             tb, nl, q :: Builder
             tb = char '\t'
             nl = char '\n'
             q  = char '?'
             
             (rndec1,rndec2) = case decodeReadName h of Just rn -> let ((y,m,d),reg,(hh,mm,ss)) = (date rn,region rn,time rn)
                                                                   in ([putDate y m d, putTime hh mm ss, putInt2 reg]
                                                                      ,[putInt (fromIntegral $ x_loc rn), putInt (fromIntegral $ y_loc rn)])
                                                        Nothing -> ([q,q,q],[q,q])
             (qleft,qright) = (clip_qual_left rh, clip_qual_right rh)
             avg_qual qs = let l = fromIntegral (L1.length qs)
                           in if l>0 then putFix 2 $ sum (map fromIntegral $ L1.unpack qs) * 100 `div` l
                             else putFix 2 0
         in mconcat $ intersperse tb ([P.fromByteString h]
                     ++ rndec1 ++ [putInt (fromIntegral qleft), putInt (fromIntegral qright)] ++ rndec2 
                     ++ [putInt (fromIntegral nb)
                        , P.fromByteString (fi $ quals $ flowgram r), P.fromByteString (fi $ quals $ flowgram tr)
                        , putInt (n_count r)
                        , avg_qual $ unQD $ quality r, avg_qual $ unQD $ quality tr]) ++ [nl]

-- ----------------------------------------------------------
-- The --filters option, summarize filters
-- ----------------------------------------------------------
sum_filters ::  [ReadBlock] -> L1.ByteString
sum_filters rs = P.toLazyByteString $ mconcat (header:map sumf1 rs)
  where 
    header = P.fromByteString $ B.pack "# name..... \tlength \tl_trim \tr_trim \tE K D M L\tSig Q20 Adp\n"
    sumf1 rb = let
      rh = read_header rb
      rn = read_name rh
      nb = fromIntegral $ num_bases rh
      (cl,cr) = (fromIntegral $ clip_qual_left rh, fromIntegral $ clip_qual_right rh)
      dfs = mconcat $ intersperse (char ' ') $
            map (\f -> if f rb then char '+' else char ' ') 
            [discard_empty, discard_key "tcag", discard_dots 0.05, discard_mixed, discard_length 186]
      tfs = mconcat $ intersperse (char ' ') $ map (\f -> putInt3 (f rb))
            [sigint, qual20 10, find_primer rapid_adapter]
      in mconcat (intersperse (char '\t') [P.fromByteString rn, putInt nb, putInt cl, putInt cr, dfs, tfs]++[char '\n'])

-- ----------------------------------------------------------
-- The -F option: Output the sequence of flows, one flow per line
-- ----------------------------------------------------------

fi :: Flow -> ByteString
fi f | f <= 9999 && f >= 0 = farray!f
     | otherwise = let (i,r) = f `divMod` 100 in B.pack (show i++"."++show r) 
     -- error ("Can't show flow values outside [0..99.99] (You had: "++show f++")")

farray :: Array Flow ByteString
farray = listArray (0,9999) [B.pack (showFFloat (Just 2) i "") | i <- [0,0.01..99.99::Double]]

tab :: ByteString
tab = B.pack "\t"

showread :: CommonHeader -> ReadBlock -> [ByteString]
showread h rd = let rh = read_header rd
                    rn = read_name rh
                    maskFlows = mask rh 1 qgroups . unpack 
                    qgroups = qgroup (B1.unpack $ flow_index rd) (map Qual $ L1.unpack $ unQD $ quality rd)
                    format p c v q = B.concat [rn,tab,B.pack (show p),tab,B.pack [c],tab,fi v,tab,B.pack (init $ drop 1 $ show q)]
                in zipWith4 format [(1::Int)..] (maskFlows $ flow h) (flowgram rd) qgroups

-- lower case based on the clip_qual values
mask :: ReadHeader -> Int -> [[a]] -> [Char] -> [Char]
mask _ _ _ [] = [] -- qgroups are infinite
mask rh p (q1:qs) (c:cs) = c' : mask rh (p+length q1) qs cs
    where c' = if fromIntegral p < clip_qual_left rh || fromIntegral p > clip_qual_right rh then toLower c else c
mask _ _ _ _ = error "internal error in 'mask'"

zipWith4 :: (a -> b -> c -> d -> e) -> [a] -> [b] -> [c] -> [d] -> [e]
zipWith4 f (a:as) (b:bs) (c:cs) (d:ds) =  f a b c d : zipWith4 f as bs cs ds
zipWith4 _ _ _ _ _ = []

-- | Take the unpacked index_offsets and quality values, and return 
--   a list of groups of quality values, each group corresponding to a flow value. 
--   Flow values < 0.5 result in empty groups.
qgroup :: [Index] -> [Qual] -> [[Qual]]
qgroup [] []       = let rest = []:rest in rest
qgroup is@(1:_) qs = let (iz,irest) = span (==0) (tail is)
                         (q1,qrest) = splitAt (length iz+1) qs
                     in q1 : qgroup irest qrest
qgroup (i:is) qs = [] : qgroup (i-1:is) qs
qgroup _ _ = error "internal error in 'qgroup'"

-- ----------------------------------------------------------
-- The -h option: Output a histogram of flow values
-- ----------------------------------------------------------

type Hist = UArray Flow Int

histogram :: String -> [[Flow]] -> [Hist]
histogram fl scores = runST $ do 
  let zero = newArray (0,9999) 0 :: ST s (STUArray s Flow Int)
  a <- zero
  c <- zero
  g <- zero
  t <- zero
  let ins1 ('A',i) = bump a i
      ins1 ('C',i) = bump c i
      ins1 ('G',i) = bump g i
      ins1 ('T',i) = bump t i
      ins1 (x,_)   = error ("Illegal character "++show x++" in flow!")
      bump ar i = readArray ar i >>= \x -> writeArray ar i (x+1)
  mapM_ ins1 (zip (cycle fl) (map (\x->if x>9999 || x<0 then 9999 else x) $ concat scores))
  a' <- unsafeFreeze a
  c' <- unsafeFreeze c
  g' <- unsafeFreeze g
  t' <- unsafeFreeze t
  return [a',c',g',t']

-- ouput a histogram in gnuplot format
showPlot :: String -> Int -> [String] -> [Hist] -> String
showPlot cmds mx hd hs = "set style data lines\nset logscale y\nset xlabel 'Flow value'\nset ylabel 'Count'\nset xtic 1\n\n"++cmds
  ++"\n\nplot "
  ++ concat (intersperse "," [" '-' ti "++show c | c <- hd]) ++ "\n"
  ++ unlines [ unlines [showFFloat (Just 2) (fromIntegral sc/100::Double) ""++"\t"++show (h!sc) | sc <- [0..fromIntegral mx]] ++ "e\n" | h <- hs]

showHist :: Int -> [String] -> [Hist] -> String
showHist mx hd hs = (if not (null hd) then concat (intersperse "\t" ("Score":hd)) ++ "\n" else "") ++
                    unlines [concat $ intersperse "\t" $ showFFloat (Just 2) (fromIntegral sc/100::Double) "" : [show (h!sc) | h <- hs] | sc <- [0..fromIntegral mx]]

histpos :: Int -> [ReadBlock] -> [Hist]
histpos mx scores = runST $ do 
  let mx' = fromIntegral mx
      zero = newArray (0,mx') 0 :: ST s (STUArray s Flow Int)
  hs <- sequence $ replicate (length $ flowgram $ head scores) zero
  let bump ar i = when (i>=0 && i<= mx') (readArray ar i >>= \x -> writeArray ar i (x+1))
  sequence_ $ concatMap (zipWith bump hs . flowgram) scores
  mapM unsafeFreeze hs
