module Print 
    (
     Builder, toLazyByteString, mconcat, fromByteString, char
    , putInt, putInt2, putInt3, putDate, putTime, putFix
    ) where 

import Data.Binary.Builder
import Data.Monoid
import Data.ByteString.Char8 (ByteString, pack)
import Data.Array.Unboxed
import Data.Char (ord)

char :: Char -> Builder
char = singleton . fromIntegral . ord

putInt :: Int -> Builder
putInt y = if y < 0 then char '-' `append` putInt' (negate y)
           else putInt' y
    where putInt' x =  case x `divMod` 1000 of
             (0,r) -> fromByteString (ints!r)
             (a,r) -> putInt a `append` putInt3 r

-- zero padded ints

putInt2 :: Int -> Builder
putInt2 x | x<100 = fromByteString (int2s!x)
          | otherwise = fromByteString (pack "xx")

putInt3 :: Int -> Builder
putInt3 x | x<1000 = fromByteString (int3s!x)
          | otherwise = fromByteString (pack "xxx")

ints, int2s, int3s :: Array Int ByteString
ints  = listArray (0,999) [pack (show i) | i <- [0..999::Int]]
int2s = listArray (0,99) (map (pack . ('0':) . show) [0..9::Int]++[pack (show i) | i <- [10..99::Int]])
int3s = listArray (0,999) (map (pack . (' ':) . (' ':) . show) [0..9::Int] ++ map (pack . (' ':) . show ) [10..99::Int] ++ map (pack . show) [100..999::Int])

putDate :: Int -> Int -> Int -> Builder
putDate y m d = mconcat [putInt y, dash, putInt2 m, dash, putInt2 d]
    where dash = char '-'

putTime :: Int -> Int -> Int -> Builder
putTime h m s = mconcat [putInt2 h, col, putInt2 m, col, putInt2 s]
    where col = char ':'

-- a bit hackish, maybe?
putFix :: Int -> Int -> Builder
putFix 1 n = let (i,f) = n `divMod` 10 in putInt i `mappend` char '.' `mappend` putInt f
putFix 2 n = let (i,f) = n `divMod` 100 in putInt i `mappend` char '.' `mappend` putInt2 f
putFix 3 n = let (i,f) = n `divMod` 1000 in putInt i `mappend` char '.' `mappend` putInt3 f
putFix _ _ = error "putFix only supports up to three fractional decimals"