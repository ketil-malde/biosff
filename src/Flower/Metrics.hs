-- Calculate various characteristics on sequence quality

module Metrics where

import Bio.Sequence.SFF
import qualified Data.ByteString.Lazy.Char8 as B

-- import Test.QuickCheck

-- | Take the fractional parts of the flows, and sum their squares (the "K²" metric)
quals :: [Flow] -> Flow
quals q = floor $ ((100::Double) - 2*(sqrt $ (/fromIntegral (length q)) $ sum $ map (fromIntegral . (^(2::Integer)) . (flip (-) 50) . (`mod` 100) . (+50)) $ q))

-- | Count number of n's in the sequence
--   The algorithm for generating Ns is a bit opaque, and appears to depend on the magnitude 
--   of the noise flow values.  We chicken out, and just count the called sequence.
n_count :: ReadBlock -> Int
n_count r = length . filter isN . clip . B.unpack . bases $ r
    where isN x = x=='N' || x == 'n'
          clip = take (right-left+1) . drop left
          right = fromIntegral $ clip_qual_right (read_header r)
          left = fromIntegral $ clip_qual_left (read_header r)
