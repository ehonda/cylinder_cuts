module Lib where

import Data.Maybe (fromMaybe)
import Numeric.Tools.Integration
import Test.LeanCheck
import Text.PrettyPrint

-- Doc Utility
----------------------------------------------------

eqDoc :: String -> String -> Doc
eqDoc name quantity = hsep $ map text
    [ name, "=", quantity ]
    
capIndent :: String -> [Doc] -> Doc
capIndent cap ds = hang (text cap) 4 $ vcat ds
    
-- Cylinder Params
----------------------------------------------------
data CylinderParams = CylinderParams
    { parL, parR, parH0, parHL :: Double }
    deriving Show
    
-- Preconditions on our parameters
--      0.0 < L < Infinity
--      0.0 < r < Infinity
--      h_0, h_L \in [0, 2r]
------------------------------------------
preCond_L :: CylinderParams -> Bool
preCond_L p =
        0.0 < len
    &&  (not $ isInfinite len)
    where
        len = parL p

preCond_r :: CylinderParams -> Bool
preCond_r p =
        0.0 < r
    &&  (not $ isInfinite r)
    where
        r = parR p

preCond_hs :: CylinderParams -> Bool
preCond_hs p =
        (0.0 <= h_0 && h_0 <= 2*r)
    &&  (0.0 <= h_L && h_L <= 2*r)
    &&  h_0 /= h_L
    where
        r = parR p
        h_0 = parH0 p
        h_L = parHL p

preCond_Params :: CylinderParams -> Bool
preCond_Params p =
        preCond_L p 
    &&  preCond_r p 
    &&  preCond_hs p
    
-- Utility functions
------------------------------------------
withHeights :: Double -> Double -> CylinderParams
withHeights h_0 h_L = CylinderParams
    { parL = 1, parR = 1, parH0 = h_0, parHL = h_L }
    
flipHeights :: CylinderParams -> CylinderParams
flipHeights p = CylinderParams 
    {   parL = parL p, parR = parR p, 
        parH0 = parHL p, parHL = parH0 p }
        
paramsDoc :: CylinderParams -> Doc
paramsDoc p = capIndent
    "Cylinder parameters:"
    [ eqDoc "L  " $ show $ parL p
    , eqDoc "r  " $ show $ parR p
    , eqDoc "h_0" $ show $ parH0 p
    , eqDoc "h_L" $ show $ parHL p
    ]
        
printParams :: CylinderParams -> IO ()
printParams p = do
    putStrLn $ render $ paramsDoc p

-- Functions
----------------------------------------------------
    
sigma :: Double -> Double -> Double
sigma r x = sqrt $ 2*r*x - x^2

heightFunction :: Double -> Double -> Double -> 
    Double -> Double
heightFunction h_0 h_L len x =
    ((len - x) * h_0 + x * h_L) / len

areaFunction :: Double -> Double -> Double
areaFunction r h = r^2 * 
    acos (1 - h / r) - (r - h) * sigma r h

volumeFormula :: CylinderParams -> Double
volumeFormula p = (r^2/c_h) * ((rho h_0) - (rho h_L) + 
    ((cr h_L)*(s h_L) - (cr h_0)*(s h_0))/3)
    where
        len = parL p
        r = parR p
        h_0 = parH0 p
        h_L = parHL p
        
        c_h = (h_L - h_0) / len
        
        rho :: Double -> Double
        rho x = (r - x) * acos((r - x)/r)
        
        cr :: Double -> Double
        cr x = 2 + ((r - x)/r)^2
        
        s :: Double -> Double
        s = sigma r

-- Comparison Formula <-> Numeric
----------------------------------------------------
numericVolume :: 
    (QuadParam -> (Double, Double) -> (Double -> Double) -> QuadRes)
    -> CylinderParams
    -> QuadRes
numericVolume f p = f defQuad (0, len) (a . h)
    where
        len = parL p
        r = parR p
        h_0 = parH0 p
        h_L = parHL p
        
        h :: Double -> Double
        h = heightFunction h_0 h_L len
        
        a :: Double -> Double
        a = areaFunction r

rombergVolume :: CylinderParams -> QuadRes
rombergVolume = numericVolume quadRomberg

trapezoidVolume :: CylinderParams -> QuadRes
trapezoidVolume = numericVolume quadTrapezoid
        
simpsonVolume :: CylinderParams -> QuadRes
simpsonVolume = numericVolume quadSimpson

resutsDoc :: CylinderParams -> Double -> QuadRes -> Doc
resutsDoc p formula estimate = vcat
    [ pDoc, fDoc, eDoc, dDoc ]
    where
        pDoc = paramsDoc p
            
        fDoc = capIndent
            "Formula result:"
            [ eqDoc "V" $ show formula ]
            
        eDoc = capIndent
            "Numeric Integration result:"
            [ eqDoc "V" $ fromMaybe ("Integration failed") $
                fmap show $ quadRes estimate
            , eqDoc "Best Estimate" $ show $
                quadBestEst estimate
            , eqDoc "Estimated Precision" $ show $
                quadPrecEst estimate
            , eqDoc "Number of Iterations" $ show $
                quadNIter estimate
            ]
            
        dDoc = capIndent
            "Difference Formula - Best Estimate:"
            [ eqDoc "V_f - V_e" $ show $ formula - 
                quadBestEst estimate
            ]

compareFormulaToNumeric :: CylinderParams -> IO ()
compareFormulaToNumeric p = do
    putStrLn $ render $ resutsDoc p formula estimate
    where
        formula = volumeFormula p
        estimate = rombergVolume p
        

-- LeanCheck
----------------------------------------------------
instance Listable CylinderParams where
    tiers = (cons4 CylinderParams) `suchThat` preCond_Params

withinTolerance :: Double -> Double -> Double -> Bool
withinTolerance delta x y = abs (x - y) < delta

prop_VolumeWithinTolerance :: Double -> 
    CylinderParams -> Bool
prop_VolumeWithinTolerance delta p =
    withinTolerance delta formula estimate
    where
        formula = volumeFormula p
        estimate = quadBestEst $ rombergVolume p

{--        
Unsuccessful LeanCheck:

Problem with the numeric integrations:

q = CylinderParams {parL = 0.4, parR = 1.5, parH0 = 0.0, parHL = 3.0}

gives NaN, but

q = CylinderParams {parL = 0.40000001, parR = 1.5, parH0 = 0.0, parHL = 3.0}

works just fine

--}
test0 = checkFor 10000 $ prop_VolumeWithinTolerance 1e-8
test1 = checkFor 100000 $ prop_VolumeWithinTolerance 1e-8