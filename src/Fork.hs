module Fork where

import Control.Concurrent
import Control.Exception

-- | Spawn a set of threads and wait for them to complete.
forkAndWait :: [IO ()] -> IO ()
forkAndWait actions = mapM myForkIO actions >>= mapM_ takeMVar
  where
    myForkIO :: IO () -> IO (MVar ())
    myForkIO io = do
      mvar <- newEmptyMVar
      _ <- forkIO (io `finally` putMVar mvar ())
      return mvar
