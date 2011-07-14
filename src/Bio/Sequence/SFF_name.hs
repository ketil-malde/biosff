module Bio.Sequence.SFF_name where

import qualified Data.ByteString.Char8 as B
import Data.ByteString.Char8 (ByteString, pack)
import Data.Array.Unboxed
import Data.Char (ord)

-- | Read names encode various information, as per this struct.
data ReadName = ReadName { date :: (Int,Int,Int)
                         , time :: (Int,Int,Int)
                         , region :: Int
                         , x_loc, y_loc :: Int } deriving Show

-- ----------------------------------------------------------
-- Decoding

decodeReadName :: ByteString -> Maybe ReadName
decodeReadName b = do t <- decodeDate $ B.take 6 b
                      r <- fst `fmap` (B.readInt $ B.take 2 $ B.drop 7 b)
                      l <- decodeLocation $ B.drop 9 b
                      return $ ReadName { date = (\[y,m,d] -> (y,m,d)) (take 3 t)
                               , time = (\[hh,mm,ss] -> (hh,mm,ss)) (drop 3 t)
                               , region = r
                               , x_loc = fst l, y_loc = snd l }

decodeLocation :: ByteString -> Maybe (Int,Int)
decodeLocation l = (`divMod` 4096) `fmap` decode36 l

decodeDate :: ByteString -> Maybe [Int]
decodeDate d    = (fixyear . reverse . (`divMods` [60,60,24,32,13])) =<< decode36 d
    where fixyear (i:is) = Just (2000+i:is)
          fixyear []     = Nothing

-- ----------------------------------------------------------
-- Encoding

encodeReadName :: ReadName -> ByteString
encodeReadName r =  B.concat [ encodeDate (date r) (time r) 
                             , encodeRegion (region r)
                             , encodeLocation (x_loc r) (y_loc r)]

encodeLocation :: Int -> Int -> ByteString
encodeLocation = undefined

encodeRegion :: Int -> ByteString
encodeRegion = undefined

encodeDate :: (Int,Int,Int) -> (Int,Int,Int) -> ByteString
encodeDate = undefined

-- ----------------------------------------------------------

divMods :: Int -> [Int] -> [Int]
divMods x (i:is) = let (a,b) = x `divMod` i
                   in b : divMods a is
divMods x [] = [x]

-- ----------------------------------------------------------
-- Decoding base36 strings

decode36 :: ByteString -> Maybe Int
decode36 s = (foldr1 (\a b -> b*36+a) . reverse) `fmap` (mapM decCh . B.unpack $ s)

{-
decode36' = dec 0
    where dec i b = case uncons b of Just (c,rest) -> dec (i*36+fromJust (decCh c)) rest
                                     Nothing       -> i
          fromJust (Just z) = z
-}

decCh :: Char -> Maybe Int
decCh x | x >= 'A' && x <= 'Z' = Just (ord x - ord 'A')
        | x >= '0' && x <= '9' = Just (26 + ord x - ord '0')
        | otherwise            = Nothing -- error ("decode36: can't decode "++show x)

encode36 :: Int -> ByteString
encode36 = pack . map (b36!) . reverse . enc
    where
      enc 0 = []
      enc i = let (a,b) = i `divMod` 36
              in b : enc a

b36 :: UArray Int Char
b36 = listArray (0,35) (['A'..'Z']++['0'..'9'])


